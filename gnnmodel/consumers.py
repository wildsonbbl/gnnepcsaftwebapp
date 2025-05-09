"consume module"

import asyncio
import base64
import io
import json
import os
import uuid
from atexit import register
from contextlib import AsyncExitStack
from typing import Any, Dict, List, Optional

import filetype
import PyPDF2
from channels.db import database_sync_to_async
from channels.generic.websocket import AsyncWebsocketConsumer
from django.conf import settings
from google.adk.agents.run_config import RunConfig
from google.adk.runners import Runner
from google.adk.sessions.in_memory_session_service import Session
from google.adk.tools.mcp_tool.mcp_toolset import (
    MCPTool,
    MCPToolset,
    StdioServerParameters,
)
from google.genai.types import Content, Part
from markdown import markdown

from . import logger
from .agents import AVAILABLE_MODELS, DEFAULT_MODEL, all_tools, gemini_models_data
from .agents_utils import get_ollama_models, is_api_key_valid, is_ollama_online
from .chat_utils import (
    APP_NAME,
    USER_ID,
    BlankLinkExtension,
    CustomJSONEncoder,
    docstring_to_html,
    session_service,
    start_agent_session,
)
from .models import ChatSession

# MCPToolset exit stack
mcp_exit_stack = AsyncExitStack()


def _cleanup_mcp_resources_on_exit():
    """Função síncrona para limpar recursos de mcp_exit_stack na saída do programa."""
    try:
        # Use asyncio.run() para executar o método assíncrono aclose().
        # Isso criará um novo loop de eventos para esta tarefa específica.
        asyncio.run(mcp_exit_stack.aclose())
        logger.info("mcp_exit_stack fechado com sucesso via atexit.")
    except RuntimeError as e:
        # Isso pode ocorrer se aclose() for chamado em uma pilha já fechada,
        # ou se houver um problema com o gerenciamento do loop de eventos na saída.
        logger.warning(
            "RuntimeError durante a limpeza atexit para mcp_exit_stack: %s. "
            "Isso pode ser normal se já estiver fechado ou devido"
            " ao estado do loop de eventos na saída.",
            e,
        )
    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.error(
            "Erro inesperado durante a limpeza atexit para mcp_exit_stack: %s", e
        )


register(_cleanup_mcp_resources_on_exit)


class CurrentChatSessionConsumer(AsyncWebsocketConsumer):
    "Current Chat Session Consumer"

    session_id: str
    runner_session: Session
    runner: Runner
    run_config = RunConfig(response_modalities=["TEXT"])
    agent_task = None
    mcp_tools: List[MCPTool] = []
    original_tools: List[Any] = all_tools
    tool_descriptions: Dict[str, str] = {
        t.__name__: docstring_to_html(t.__doc__) or "No description available"
        for t in all_tools
    }


class CurrentChatSessionConsumerUtils(CurrentChatSessionConsumer):
    "Chat Consumer Utils"

    async def load_session_data(
        self,
        session: ChatSession,
        current_tool_map: Dict[str, Any],
        valid_selected_tools: List[str],
    ):
        "load session data"
        # Initialize the agent with the session
        self.runner, self.runner_session = await start_agent_session(
            self.session_id,
            session.model_name,
            tools=[current_tool_map[tool] for tool in valid_selected_tools],
        )

        # Send the current session info to the client
        await self.send(
            text_data=json.dumps(
                {
                    "action": "session_loaded",
                    "session_id": self.session_id,
                    "name": session.name,
                    "model_name": session.model_name,
                    "available_models": AVAILABLE_MODELS,
                    "available_tools": list(current_tool_map),
                    "selected_tools": valid_selected_tools,
                    "tool_descriptions": self.tool_descriptions,
                    "mcp_config_path": str(settings.MCP_SERVER_CONFIG),
                },
                cls=CustomJSONEncoder,
            )
        )

        await self.send(
            text_data=json.dumps(
                {"action": "load_messages", "messages": session.messages},
                cls=CustomJSONEncoder,
            )
        )

    async def validate_and_update_tools(
        self, session: ChatSession, current_tool_map: Dict[str, Any]
    ) -> List[str]:
        "validate and update tools"
        valid_selected_tools = [
            tool for tool in session.selected_tools if tool in current_tool_map
        ]
        if len(valid_selected_tools) != len(session.selected_tools):
            # Update the session if some selected tools are no
            # longer valid (e.g., after MCP activation)
            await database_sync_to_async(
                ChatSession.objects.filter(session_id=self.session_id).update
            )(selected_tools=valid_selected_tools)
        return valid_selected_tools

    def get_current_tool_map(self, current_tools: List[Any]) -> Dict[str, Any]:
        """get current tool map"""
        current_tool_map = {}
        for t in current_tools:
            if hasattr(t, "__name__"):  # Check for regular tools
                current_tool_map[t.__name__] = t
            elif hasattr(t, "name"):  # Check for MCP tools
                current_tool_map[t.name] = t
        return current_tool_map

    async def send_error_message(self, error_text):
        """Sends an error message formatted as a system message."""
        error_message = {
            "msg": markdown(
                f"**Error:** {error_text}", extensions=[BlankLinkExtension()]
            ),
            "source": "assistant",  # Or a dedicated 'system' source
        }
        await self.send(text_data=json.dumps({"text": error_message}))

    async def get_mcp(
        self, command: str, args: List[str], env: Optional[Dict[str, str]] = None
    ):
        "get mcp toolset"

        tools, _ = await MCPToolset.from_server(
            connection_params=StdioServerParameters(
                command=command, args=args, env=env
            ),
            async_exit_stack=mcp_exit_stack,
        )

        return tools

    async def activate_mcp_server(self):
        "activate the servers on the config"
        # Clear existing MCP tools before activation
        self.mcp_tools = []
        activated_tool_names = []
        error_message = None

        try:
            with open(settings.MCP_SERVER_CONFIG, "r", encoding="utf-8") as config_file:
                mcp_server_config: Dict[str, Dict[str, Dict[str, Any]]] = json.load(
                    config_file
                )
        except FileNotFoundError:
            error_message = "MCP configuration file not found, config file first."
            logger.error(error_message)
            mcp_server_config = {}
        except json.JSONDecodeError as e:
            error_message = (
                f"Error decoding JSON from MCP configuration"
                f" file: {settings.MCP_SERVER_CONFIG}. Error: {e}"
            )
            logger.error(error_message)
            mcp_server_config = {}
        except (  # pylint: disable=broad-exception-caught
            Exception
        ) as e:  # Catch other potential errors during file opening/reading
            error_message = (
                f"Error reading MCP configuration"
                f" file: {settings.MCP_SERVER_CONFIG}. Error: {e}"
            )
            logger.error(error_message)
            mcp_server_config = {}

        if not error_message and "mcpServers" in mcp_server_config:
            logger.debug(mcp_server_config)
            for mcpserver_name in mcp_server_config["mcpServers"]:
                mcpserver = mcp_server_config["mcpServers"][mcpserver_name]
                command = mcpserver.get("command")
                args = mcpserver.get("args")
                env = mcpserver.get("env")

                # Basic validation for command and args
                if not isinstance(command, str) or not command:
                    logger.error(
                        "Invalid or missing 'command' for MCP server: %s",
                        mcpserver_name,
                    )
                    continue  # Skip this server
                if not isinstance(args, list):
                    logger.error(
                        "Invalid or missing 'args' for MCP server: %s. Must be a list.",
                        mcpserver_name,
                    )
                    continue  # Skip this server

                try:
                    new_tools = await self.get_mcp(command, args, env)
                    self.mcp_tools.extend(new_tools)
                    activated_tool_names.extend([t.name for t in new_tools])
                    for t in new_tools:
                        self.tool_descriptions[t.name] = t.description
                    logger.info(
                        "Activated MCP tools from server '%s': %s",
                        mcpserver_name,
                        [t.name for t in new_tools],
                    )
                except Exception as e:  # pylint: disable=broad-exception-caught
                    error_message = (
                        f"Failed to activate MCP server '{mcpserver_name}': {e}"
                    )
                    logger.error(error_message)
                    # Decide if we should stop all activation on first error, or continue
                    # For now, let's record the error and continue
                    # break # Uncomment to stop on first error

        return activated_tool_names, error_message

    async def process_uploaded_file(self, file_info: dict):
        """Processes base64 encoded file data into a Gemini Part."""
        try:
            # Decode base64 data URL
            # Format is "data:[<mediatype>][;base64],<data>"
            header, encoded_data = file_info["data"].split(",", 1)
            # data:image/png;base64 -> image/png
            mime_type = header.split(":")[1].split(";")[0]
            file_bytes = base64.b64decode(encoded_data)
            file_name = file_info.get("name", "uploaded_file")
            kind = filetype.guess(file_bytes)
            if kind is None:
                logger.warning(
                    "Could not determine file type for '%s' using filetype.",
                    file_name,
                )
                raise ValueError(
                    f"Unsupported or unrecognized file type for '{file_name}'."
                )

            if kind.mime != mime_type:
                logger.warning(
                    "Client provided MIME type '%s' but filetype detected '%s' for file '%s'. "
                    "Using detected type.",
                    mime_type,
                    kind.mime,
                    file_name,
                )
                mime_type = kind.mime

            processed_file_info_for_db = {
                "name": file_name,
                "type": mime_type,
                # Store original data URL for potential image display
                "data": file_info["data"] if mime_type.startswith("image/") else None,
            }

            if mime_type == "application/pdf":
                # Extract text from PDF
                pdf_reader = PyPDF2.PdfReader(io.BytesIO(file_bytes))
                pdf_text = f"Content from PDF '{file_name}':\n"
                for page in pdf_reader.pages:
                    pdf_text += page.extract_text() + "\n"
                # Create a text Part for PDF content
                return Part.from_text(text=pdf_text), processed_file_info_for_db
            if mime_type.startswith("image/"):
                # Create an image Part
                return (
                    Part.from_bytes(data=file_bytes, mime_type=mime_type),
                    processed_file_info_for_db,
                )

            # Handle other file types (e.g., treat as plain text or reject)
            # For simplicity, let's try decoding as text, warning it might fail
            try:
                text_content = file_bytes.decode("utf-8")
                logger.warning("Treating file '%s' (%s) as text.", file_name, mime_type)
                return (
                    Part.from_text(
                        text=f"Content from file '{file_name}':\n{text_content}"
                    ),
                    processed_file_info_for_db,
                )
            except UnicodeDecodeError as exc:
                raise ValueError(
                    f"Unsupported file type: {mime_type}. Could not decode as text."
                ) from exc

        except Exception as e:
            logger.error(
                "Error decoding/processing file %s: %s", file_info.get("name"), e
            )
            raise  # Re-raise the exception to be caught in receive()


class ChatSessionsDBOperations(CurrentChatSessionConsumerUtils):
    "Session DB handling"

    async def get_or_create_last_session(self) -> ChatSession:
        "get or create the last session"
        last_session = await self.get_last_session()
        if last_session:
            self.session_id = str(last_session.session_id)
            session = last_session
        else:
            # Create a new session
            self.session_id = str(uuid.uuid4())
            session = await self.get_or_create_session(self.session_id)
        return session

    @database_sync_to_async
    def get_last_session(self) -> Optional[ChatSession]:
        """Get the most recently updated session"""
        try:
            return ChatSession.objects.order_by("-updated_at").first()
        except ChatSession.DoesNotExist:
            return None

    @database_sync_to_async
    def get_or_create_session(
        self,
        session_id,
        name="New Session",
        model_name=DEFAULT_MODEL,
        selected_tools=None,
    ) -> ChatSession:
        """Get or create a session"""
        try:
            # Try to get existing session
            return ChatSession.objects.get(session_id=session_id)
        except ChatSession.DoesNotExist:
            # Create new session
            if selected_tools is None:
                selected_tools = [t.__name__ for t in all_tools]
            return ChatSession.objects.create(
                session_id=session_id,
                name=name,
                model_name=model_name,
                selected_tools=selected_tools,
            )

    @database_sync_to_async
    def save_message_to_db(self, message):
        """Save message to database"""
        session = ChatSession.objects.get(session_id=self.session_id)
        session.add_message(message)
        # No need to return session here unless used immediately after

    def serialize_session(self, session):
        """Convert session data to JSON-serializable format"""
        return {
            "session_id": str(session["session_id"]),  # Convert UUID to string
            "name": session["name"],
            "created_at": session[
                "created_at"
            ].isoformat(),  # Convert datetime to ISO format string
            "updated_at": session["updated_at"].isoformat(),
        }

    @database_sync_to_async
    def get_all_sessions(self):
        """Get all sessions"""
        sessions = list(
            ChatSession.objects.values("session_id", "name", "created_at", "updated_at")
        )
        # Convert each session to a JSON-serializable format
        return [self.serialize_session(session) for session in sessions]

    @database_sync_to_async
    def delete_session(self, session_id):
        """Delete a session"""
        try:
            session = ChatSession.objects.get(session_id=session_id)
            session.delete()

            session_service.delete_session(
                app_name=APP_NAME, user_id=USER_ID, session_id=session_id
            )

            return True
        except ChatSession.DoesNotExist:
            return False


class ChatConsumerMessagingOperations(ChatSessionsDBOperations):
    "messaging functions"

    async def agent_to_client_messaging(self, text, file_part=None):
        """Agent to client communication"""
        try:
            if not await self.handle_api_key_validation():
                return
            parts = [Part.from_text(text=text)]
            if file_part:
                parts.append(file_part)

            await self.process_user_message(parts)

        except asyncio.CancelledError:
            # Task foi cancelada (usuário clicou em Parar)
            pass
        except Exception as e:  # pylint: disable=broad-exception-caught
            logger.error(e)
            message = {
                "msg": markdown(
                    f"***Error with agent**: `{e}`*", extensions=[BlankLinkExtension()]
                ),
                "source": "assistant",
            }
            await self.send(text_data=json.dumps({"text": message}))
            await self.save_message_to_db(message)
        finally:
            await self.send(
                text_data=json.dumps({"action": "end_turn", "type": "end_of_turn"})
            )

    async def handle_api_key_validation(self) -> bool:
        "handle api key validation"
        current_session: ChatSession = await self.get_or_create_session(self.session_id)
        model_name = current_session.model_name
        key_is_valid = True

        # Check for API key if using a Gemini model
        if model_name and model_name.lower().startswith("gemini"):
            google_api_key = os.getenv("GOOGLE_API_KEY", os.getenv("GEMINI_API_KEY"))
            key_is_valid = False
            if google_api_key:
                # is_api_key_valid can be slow as it makes a network request
                # It's called here to ensure the key is actually functional
                key_is_valid = is_api_key_valid(google_api_key)

            if not key_is_valid:
                await self.send(
                    text_data=json.dumps(
                        {
                            "action": "api_key_required",
                            "model_name": model_name,
                        }
                    )
                )
                error_message_for_log = {
                    "msg": markdown(
                        f"**System Message:** Cannot process message. "
                        f"Google API Key is missing or invalid for model `{model_name}`. "
                        "Please ensure the `GOOGLE_API_KEY` or `GEMINI_API_KEY` "
                        "environment variable is correctly set. "
                        "You can find more about Gemini API Keys "
                        "<a href='https://ai.google.dev/gemini-api/docs/api-key' "
                        "target='_blank'>here</a>.",
                        extensions=[BlankLinkExtension()],
                    ),
                    "source": "assistant",
                }
                await self.send(text_data=json.dumps({"text": error_message_for_log}))
                await self.save_message_to_db(error_message_for_log)
                await self.send(
                    text_data=json.dumps(
                        {"action": "end_turn", "type": "api_key_error"}
                    )
                )
        return key_is_valid

    async def client_to_agent_messaging(self, text, file_info_for_db=None):
        """Client to agent communication"""
        # Store user message in database
        user_message = {
            "msg": markdown(text, extensions=[BlankLinkExtension()]),
            "source": "user",
        }
        if file_info_for_db:
            # Add file info for display in chat log
            user_message["file_info"] = file_info_for_db

        await self.save_message_to_db(user_message)
        await self.send(
            text_data=json.dumps({"text": user_message}),
        )
        await self.send(text_data=json.dumps({"action": "ongoing_turn"}))

    async def bot_to_client_messaging(self):
        """Bot to client communication, for testing"""
        textes = [
            "***Error with agent**: `{error}`*",
            "booot 2",
            "booot 3",
            "'I'M A BOT",
        ]
        for text in textes:
            await self.send(
                text_data=json.dumps(
                    {
                        "text": {
                            "msg": markdown(text, extensions=[BlankLinkExtension()]),
                            "source": "assistant",
                        }
                    }
                ),
            )
            await asyncio.sleep(2)
        await self.send(
            text_data=json.dumps(
                {
                    "text": {
                        "msg": "Turn completed, brother",
                        "source": "turn_end",
                    }
                }
            ),
        )

    async def process_user_message(self, parts):
        "process user message"
        async for event in self.runner.run_async(
            new_message=Content(role="user", parts=parts),
            user_id=USER_ID,
            session_id=self.session_id,
            run_config=self.run_config,
        ):
            if event.interrupted:
                await self.send(
                    text_data=json.dumps({"action": "end_turn", "type": "interrupted"})
                )
                break

            all_parts = event.content and event.content.parts
            if not all_parts or event.partial:
                continue
            all_texts = ""
            for part in all_parts:
                if part.text:
                    resp_text = part.text  # Renamed to avoid conflict
                    all_texts += resp_text + "\n\n"
                elif part.function_call and part.function_call.name:
                    resp_text = "**Calling:** `" + part.function_call.name + "`"
                    all_texts += resp_text + "\n\n"
                else:
                    resp_text = None
                if resp_text:
                    message = {
                        "msg": markdown(resp_text, extensions=[BlankLinkExtension()]),
                        "source": "assistant",
                    }
                    await self.send(text_data=json.dumps({"text": message}))
                    await self.send(text_data=json.dumps({"action": "ongoing_turn"}))
                    await self.save_message_to_db(message)

                else:
                    all_texts += f"No text or function_call in part: {part}"
            await asyncio.sleep(0.5)


class ChatConsumerHandleActions(ChatConsumerMessagingOperations):
    "handle actions class"

    async def handle_actions(self, text_data_json: Dict[str, str]):
        """Handle actions such as creating a new session or deleting a session"""
        action = text_data_json["action"]
        action_handlers = {
            "change_model": self.handle_change_model,
            "delete_session": self.handle_delete_session,
            "get_sessions": self.handle_get_sessions,
            "create_session": self.handle_create_session,
            "load_session": self.handle_load_session,
            "rename_session": self.handle_rename_session,
            "stop_generating": self.handle_stop_generating,
            "change_tools": self.handle_change_tools,
            "activate_mcp": self.handle_activate_mcp,
            "get_mcp_config_content": self.handle_get_mcp_config_content,
            "save_mcp_config_content": self.handle_save_mcp_config_content,
            "update_ollama_models": self.handle_update_ollama_models,
        }

        handler = action_handlers.get(action)
        if handler:
            # For handlers that don't take text_data_json
            if action in [
                "get_sessions",
                "stop_generating",
                "activate_mcp",
                "get_mcp_config_content",
                "update_ollama_models",
            ]:
                await handler()
            else:
                await handler(text_data_json)
        else:
            logger.warning("Unknown action received: %s", action)
            # Optionally, send an error message back to the client
            # await self.send_error_message(f"Unknown action: {action}")

    async def handle_update_ollama_models(self):
        """Checks if Ollama is online and updates the list of available models."""
        if not is_ollama_online():
            await self.send(text_data=json.dumps({"action": "ollama_offline"}))
            return

        # Start with the initial static list of models (primarily Gemini)
        # This is the static list defined in agents.py before dynamic loading

        refreshed_available_models = AVAILABLE_MODELS
        if (
            gemini_models_data
            and "models" in gemini_models_data
            and gemini_models_data["models"]
        ):
            # Override with fetched Gemini models if available and not empty
            refreshed_available_models = gemini_models_data["models"]

        # Now, add Ollama models
        ollama_models_data = get_ollama_models()
        if (
            ollama_models_data
            and "models" in ollama_models_data
            and ollama_models_data.get("models")
        ):  # Check if models list is not empty
            ollama_model_names = [
                f"ollama_chat/{model['name']}" for model in ollama_models_data["models"]
            ]
            refreshed_available_models.extend(ollama_model_names)

        await self.send(
            text_data=json.dumps(
                {
                    "action": "available_models_updated",
                    "available_models": refreshed_available_models,
                }
            )
        )

    async def handle_activate_mcp(self):
        "handle activate mcp action"
        activated_tool_names, error_message = await self.activate_mcp_server()

        if error_message:
            await self.send(
                text_data=json.dumps(
                    {
                        "action": "mcp_activation_failed",
                        "error": error_message,
                    }
                )
            )
        else:
            # Successfully activated (or no servers defined), re-initialize agent with new tools
            session: ChatSession = await self.get_or_create_session(self.session_id)
            current_tools = self.original_tools + self.mcp_tools
            current_tool_map = self.get_current_tool_map(current_tools)

            valid_selected_tools = await self.validate_and_update_tools(
                session, current_tool_map
            )

            self.runner, self.runner_session = await start_agent_session(
                self.session_id,
                session.model_name,
                tools=[current_tool_map[t] for t in valid_selected_tools],
            )

            await self.send(
                text_data=json.dumps(
                    {
                        "action": "mcp_activated",
                        "activated_tools": activated_tool_names,
                        "available_tools": list(current_tool_map),
                        "selected_tools": valid_selected_tools,
                        "tool_descriptions": self.tool_descriptions,
                    }
                )
            )

    async def handle_change_tools(self, text_data_json):
        "handle change tools"
        current_tools = self.original_tools + self.mcp_tools
        current_tool_map = self.get_current_tool_map(current_tools)
        valid_selected_tools = [
            t for t in text_data_json["tools"] if t in current_tool_map
        ]

        # update db
        await database_sync_to_async(
            ChatSession.objects.filter(session_id=self.session_id).update
        )(selected_tools=valid_selected_tools)

        session: ChatSession = await self.get_or_create_session(self.session_id)
        self.runner, self.runner_session = await start_agent_session(
            self.session_id,
            session.model_name,
            tools=[current_tool_map[t] for t in session.selected_tools],
        )
        await self.send(
            text_data=json.dumps(
                {
                    "action": "tools_changed",
                    "selected_tools": session.selected_tools,
                    "available_tools": list(current_tool_map),
                }
            )
        )

    async def handle_stop_generating(self):
        "handle stop generating"
        if self.agent_task and not self.agent_task.done():
            self.agent_task.cancel()
        await self.send(
            text_data=json.dumps({"action": "end_turn", "type": "stop_action"})
        )

    async def handle_rename_session(self, text_data_json):
        "handle rename session"
        session_id = text_data_json["session_id"]
        name = text_data_json["name"]
        await database_sync_to_async(
            ChatSession.objects.filter(session_id=session_id).update
        )(name=name)

        await self.send(
            text_data=json.dumps(
                {
                    "action": "session_renamed",
                    "session_id": session_id,
                    "name": name,
                }
            )
        )

    async def handle_load_session(self, text_data_json):
        "handle load session"
        self.session_id = text_data_json["session_id"]

        session = await self.get_or_create_session(self.session_id)
        current_tools = self.original_tools + self.mcp_tools
        current_tool_map = self.get_current_tool_map(current_tools)
        valid_selected_tools = await self.validate_and_update_tools(
            session, current_tool_map
        )

        await self.load_session_data(session, current_tool_map, valid_selected_tools)

    async def handle_create_session(self, text_data_json: Dict[str, str]):
        "handle create session"
        self.session_id = str(uuid.uuid4())
        name = text_data_json.get("name", "New Session")
        model_name = text_data_json.get("model_name", DEFAULT_MODEL)
        current_tools = self.original_tools + self.mcp_tools
        current_tool_map = self.get_current_tool_map(current_tools)
        selected_tools_names = text_data_json.get(
            "tools",
            list(current_tool_map),  # Default to all current tools
        )
        valid_selected_tools = [
            name for name in selected_tools_names if name in current_tool_map
        ]

        session = await self.get_or_create_session(
            self.session_id,
            name=name,
            model_name=model_name,
            selected_tools=valid_selected_tools,
        )

        await self.send(
            text_data=json.dumps(
                {
                    "action": "session_created",
                    "session_id": self.session_id,
                    "name": session.name,
                }
            )
        )

        await self.load_session_data(session, current_tool_map, valid_selected_tools)

    async def handle_get_sessions(self):
        "handle get sessions"
        sessions = await self.get_all_sessions()

        await self.send(
            text_data=json.dumps({"action": "sessions_list", "sessions": sessions})
        )

    async def handle_delete_session(self, text_data_json):
        "handle delete session"
        session_id = text_data_json["session_id"]
        success = await self.delete_session(session_id)

        # If the deleted session was the current one, load the most recent session
        if success and session_id == self.session_id:
            session = await self.get_or_create_last_session()

            current_tools = self.original_tools + self.mcp_tools
            current_tool_map = self.get_current_tool_map(current_tools)
            valid_selected_tools = await self.validate_and_update_tools(
                session, current_tool_map
            )
            await self.load_session_data(
                session, current_tool_map, valid_selected_tools
            )

        # Send confirmation and refresh sessions list
        await self.send(
            text_data=json.dumps(
                {
                    "action": "session_deleted",
                    "success": success,
                    "session_id": session_id,
                }
            )
        )

        # Send updated sessions list
        sessions = await self.get_all_sessions()
        await self.send(
            text_data=json.dumps(
                {"action": "sessions_list", "sessions": sessions},
                cls=CustomJSONEncoder,
            )
        )

    async def handle_change_model(self, text_data_json):
        "handle change model"
        model_name = text_data_json["model_name"]
        if model_name in AVAILABLE_MODELS:
            # Update the session with the new model
            await database_sync_to_async(
                ChatSession.objects.filter(session_id=self.session_id).update
            )(model_name=model_name)

            current_tools = self.original_tools + self.mcp_tools
            current_tool_map = self.get_current_tool_map(current_tools)
            valid_selected_tools = [
                name for name in text_data_json["tools"] if name in current_tool_map
            ]

            self.runner, self.runner_session = await start_agent_session(
                self.session_id,
                model_name,
                [current_tool_map[name] for name in valid_selected_tools],
            )

            await self.send(
                text_data=json.dumps(
                    {
                        "action": "model_changed",
                        "model_name": model_name,
                    }
                )
            )

            # Add a system message to indicate model change
            system_message = {
                "msg": markdown(
                    f"Model changed to `{model_name}`",
                    extensions=[BlankLinkExtension()],
                ),
                "source": "user",
            }
            await self.save_message_to_db(system_message)
            await self.send(text_data=json.dumps({"text": system_message}))

    async def handle_text(self, text, file_info):
        "Handle received text"
        file_part = None
        processed_file_info_for_db = None  # Store info for DB/frontend display
        if file_info:
            try:
                file_part, processed_file_info_for_db = (
                    await self.process_uploaded_file(file_info)
                )
            except Exception as e:  # pylint: disable=broad-exception-caught
                logger.error("Error processing file: %s", e)
                # Send error message back to client? (Optional)
                await self.send_error_message(
                    f"Failed to process file: {file_info.get('name', 'Unknown')}. Error: {e}"
                )
                return  # Stop processing this message

        await self.client_to_agent_messaging(text, processed_file_info_for_db)
        if self.agent_task and not self.agent_task.done():
            self.agent_task.cancel()
        self.agent_task = asyncio.create_task(
            self.agent_to_client_messaging(text, file_part)
        )

    async def handle_get_mcp_config_content(self):
        """Handles request to get MCP server configuration content."""
        try:
            with open(settings.MCP_SERVER_CONFIG, "r", encoding="utf-8") as f:
                content = f.read()
            await self.send(
                text_data=json.dumps(
                    {"action": "mcp_config_content", "content": content}
                )
            )
        except FileNotFoundError:
            logger.error(
                "MCP configuration file not found at: %s", settings.MCP_SERVER_CONFIG
            )
            await self.send(
                text_data=json.dumps(
                    {
                        "action": "mcp_config_content",
                        "error": (
                            "MCP configuration file not found."
                            " The file will be created on save."
                        ),
                    }
                )
            )
        except Exception as e:  # pylint: disable=broad-exception-caught
            logger.error("Error reading MCP configuration file: %s", e)
            await self.send(
                text_data=json.dumps({"action": "mcp_config_content", "error": str(e)})
            )

    async def handle_save_mcp_config_content(self, text_data_json: Dict[str, str]):
        """Handles request to save MCP server configuration content."""
        new_content = text_data_json.get("content")
        if new_content is None:
            await self.send(
                text_data=json.dumps(
                    {
                        "action": "mcp_config_saved",
                        "success": False,
                        "error": "No content provided.",
                    }
                )
            )
            return

        try:
            # Validate if the content is valid JSON before writing
            json.loads(new_content)
            with open(settings.MCP_SERVER_CONFIG, "w", encoding="utf-8") as f:
                f.write(new_content)
            logger.info(
                "MCP configuration file updated: %s", settings.MCP_SERVER_CONFIG
            )
            await self.send(
                text_data=json.dumps({"action": "mcp_config_saved", "success": True})
            )
        except json.JSONDecodeError as e:
            logger.error("Error decoding new MCP config content as JSON: %s", e)
            await self.send(
                json.dumps(
                    {
                        "action": "mcp_config_saved",
                        "success": False,
                        "error": f"Invalid JSON format: {e}",
                    }
                )
            )
        except Exception as e:  # pylint: disable=broad-exception-caught
            logger.error("Error writing MCP configuration file: %s", e)
            await self.send(
                json.dumps(
                    {
                        "action": "mcp_config_saved",
                        "success": False,
                        "error": str(e),
                    }
                )
            )


class ChatConsumer(ChatConsumerHandleActions):
    "Chat consumer"

    async def connect(self):
        """Connect to the websocket"""
        await self.accept()

        # Get the last session or create a new one if none exists
        session = await self.get_or_create_last_session()

        # Combine original tools and any already activated MCP tools
        current_tools = self.original_tools + self.mcp_tools
        current_tool_map = self.get_current_tool_map(current_tools)

        # Filter selected tools based on currently available tools
        valid_selected_tools = await self.validate_and_update_tools(
            session, current_tool_map
        )

        await self.load_session_data(session, current_tool_map, valid_selected_tools)

    async def disconnect(self, code):
        await mcp_exit_stack.aclose()

    async def receive(self, text_data=None, bytes_data=None):
        assert text_data is not None
        text_data_json: Dict[str, str] = json.loads(text_data)

        # Handle different types of messages
        if "action" in text_data_json:
            await self.handle_actions(text_data_json)

        elif "text" in text_data_json:
            text = text_data_json["text"]
            file_info = text_data_json.get("file")  # Get potential file info

            # Basic validation
            if not text.strip() and not file_info:
                return  # Ignore empty messages without files

            # Process file if present
            await self.handle_text(text, file_info)
