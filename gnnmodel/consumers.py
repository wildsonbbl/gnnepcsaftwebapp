"consume module"

import asyncio
import base64
import io
import json
import uuid
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
from google.adk.tools.mcp_tool.mcp_toolset import MCPToolset, StdioServerParameters
from google.genai.types import Content, Part
from markdown import markdown

from . import logger
from .agents import AVAILABLE_MODELS, DEFAULT_MODEL, all_tools
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

tool_map = {t.__name__: t for t in all_tools}
tool_descriptions = {
    t.__name__: docstring_to_html(t.__doc__) or "No description available"
    for t in all_tools
}

# MCPToolset exit stack
mcp_exit_stack = AsyncExitStack()


class ChatConsumer(AsyncWebsocketConsumer):
    "Chat consumer"

    session_id: str
    runner_session: Session
    runner: Runner
    run_config = RunConfig(response_modalities=["TEXT"])
    agent_task = None
    mcp_tools: List = []

    async def connect(self):
        """Connect to the websocket"""
        await self.accept()

        # Get the last session or create a new one if none exists
        last_session = await self.get_last_session()
        if last_session:
            self.session_id = str(last_session.session_id)
            session = last_session
        else:
            # Create a new session
            self.session_id = str(uuid.uuid4())
            session = await self.get_or_create_session(self.session_id)
        assert isinstance(session, ChatSession)
        # Initialize the agent with the session
        self.runner, self.runner_session = await start_agent_session(
            self.session_id,
            session.model_name,
            tools=[tool_map[tool] for tool in session.selected_tools] + self.mcp_tools,
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
                    "available_tools": [t.__name__ for t in all_tools],
                    "selected_tools": session.selected_tools,
                    "tool_descriptions": tool_descriptions,
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

    async def disconnect(self, code):
        """
        Called when a WebSocket connection is closed.
        """
        await mcp_exit_stack.aclose()

    async def receive(self, text_data=None, bytes_data=None):
        assert text_data is not None
        text_data_json = json.loads(text_data)
        assert isinstance(text_data_json, dict)

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

        # Load configuration from the file path specified in settings

        try:
            with open(settings.MCP_SERVER_CONFIG, "r", encoding="utf-8") as config_file:
                mcp_server_config: Dict[str, Dict[str, Dict[str, Any]]] = json.load(
                    config_file
                )
        except FileNotFoundError:
            logger.error(
                "MCP configuration file not found at: %s", settings.MCP_SERVER_CONFIG
            )
            mcp_server_config = {}
        except json.JSONDecodeError:
            logger.error(
                "Error decoding JSON from MCP configuration file: %s",
                settings.MCP_SERVER_CONFIG,
            )
            mcp_server_config = {}

        logger.debug(mcp_server_config)

        if "mcpServers" in mcp_server_config:
            for mcpserver_name in mcp_server_config["mcpServers"]:
                mcpserver = mcp_server_config["mcpServers"][mcpserver_name]
                command = mcpserver.get("command")
                args = mcpserver.get("args")
                assert isinstance(command, str)
                assert isinstance(args, List)
                env = mcpserver.get("env")
                self.mcp_tools += await self.get_mcp(command, args, env)

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

    @database_sync_to_async
    def get_last_session(self):
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
    ):
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

    async def handle_actions(self, text_data_json: dict):  # pylint: disable=R0915,R0912
        """Handle actions such as creating a new session or deleting a session"""
        action = text_data_json["action"]
        if action == "change_model":
            model_name = text_data_json["model_name"]
            if model_name in AVAILABLE_MODELS:
                # Update the session with the new model
                await database_sync_to_async(
                    ChatSession.objects.filter(session_id=self.session_id).update
                )(model_name=model_name)

                selected_tools = [
                    tool_map[name]
                    for name in text_data_json["tools"]
                    if name in tool_map
                ]
                self.runner, self.runner_session = await start_agent_session(
                    self.session_id, model_name, selected_tools + self.mcp_tools
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

        elif action == "delete_session":
            session_id = text_data_json["session_id"]
            success = await self.delete_session(session_id)

            # If the deleted session was the current one, load the most recent session
            if success and session_id == self.session_id:
                last_session = await self.get_last_session()
                if last_session:
                    self.session_id = str(last_session.session_id)
                    session = last_session
                else:
                    # Create a new session
                    self.session_id = str(uuid.uuid4())
                    session = await self.get_or_create_session(self.session_id)
                assert isinstance(session, ChatSession)

                self.runner, self.runner_session = await start_agent_session(
                    self.session_id,
                    session.model_name,
                    tools=[tool_map[name] for name in session.selected_tools]
                    + self.mcp_tools,
                )

                await self.send(
                    text_data=json.dumps(
                        {
                            "action": "session_loaded",
                            "session_id": self.session_id,
                            "name": session.name,
                            "model_name": session.model_name,
                            "available_models": AVAILABLE_MODELS,
                            "available_tools": [t.__name__ for t in all_tools],
                            "selected_tools": session.selected_tools,
                        },
                        cls=CustomJSONEncoder,
                    )
                )

                await self.send(
                    text_data=json.dumps(
                        {"action": "load_messages", "messages": session.messages}
                    )
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

        elif action == "get_sessions":
            # Return list of available sessions
            sessions = await self.get_all_sessions()

            await self.send(
                text_data=json.dumps({"action": "sessions_list", "sessions": sessions})
            )

        elif action == "create_session":
            # Create a new session
            self.session_id = str(uuid.uuid4())
            name = text_data_json.get("name", "New Session")
            model_name = text_data_json.get("model_name", DEFAULT_MODEL)
            selected_tools = text_data_json.get(
                "tools", [t.__name__ for t in all_tools]
            )
            session = await self.get_or_create_session(
                self.session_id,
                name=name,
                model_name=model_name,
                selected_tools=selected_tools,
            )
            assert isinstance(session, ChatSession)
            self.runner, self.runner_session = await start_agent_session(
                self.session_id,
                model_name=session.model_name,
                tools=[tool_map[t] for t in session.selected_tools] + self.mcp_tools,
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

            await self.send(
                text_data=json.dumps(
                    {
                        "action": "session_loaded",
                        "session_id": self.session_id,
                        "name": session.name,
                        "model_name": session.model_name,
                        "available_models": AVAILABLE_MODELS,
                        "available_tools": [t.__name__ for t in all_tools],
                        "selected_tools": session.selected_tools,
                    },
                    cls=CustomJSONEncoder,
                )
            )

        elif action == "load_session":
            # Load an existing session
            self.session_id = text_data_json["session_id"]

            session = await self.get_or_create_session(self.session_id)
            assert isinstance(session, ChatSession)
            self.runner, self.runner_session = await start_agent_session(
                self.session_id,
                session.model_name,
                tools=[tool_map[t] for t in session.selected_tools] + self.mcp_tools,
            )

            await self.send(
                text_data=json.dumps(
                    {
                        "action": "session_loaded",
                        "session_id": self.session_id,
                        "name": session.name,
                        "model_name": session.model_name,
                        "available_models": AVAILABLE_MODELS,
                        "available_tools": [t.__name__ for t in all_tools],
                        "selected_tools": session.selected_tools,
                    },
                    cls=CustomJSONEncoder,
                )
            )

            await self.send(
                text_data=json.dumps(
                    {"action": "load_messages", "messages": session.messages}
                )
            )

        elif action == "rename_session":
            # Rename a session
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

        elif action == "stop_generating":
            if self.agent_task and not self.agent_task.done():
                self.agent_task.cancel()
            await self.send(
                text_data=json.dumps({"action": "end_turn", "type": "stop_action"})
            )

        elif action == "change_tools":
            selected_tools = [
                tool_map[t] for t in text_data_json["tools"] if t in tool_map
            ]
            await database_sync_to_async(
                ChatSession.objects.filter(session_id=self.session_id).update
            )(selected_tools=[t.__name__ for t in selected_tools])

            session = await self.get_or_create_session(self.session_id)
            self.runner, self.runner_session = await start_agent_session(
                self.session_id,
                session.model_name,
                tools=selected_tools + self.mcp_tools,
            )
            await self.send(
                text_data=json.dumps(
                    {
                        "action": "tools_changed",
                        "selected_tools": [t.__name__ for t in selected_tools],
                    }
                )
            )

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

    async def agent_to_client_messaging(self, text, file_part=None):
        """Agent to client communication"""
        try:
            parts = [Part.from_text(text=text)]
            if file_part:
                parts.append(file_part)

            async for event in self.runner.run_async(
                new_message=Content(role="user", parts=parts),
                user_id=USER_ID,
                session_id=self.session_id,
                run_config=self.run_config,
            ):

                if event.interrupted:
                    await self.send(
                        text_data=json.dumps(
                            {"action": "end_turn", "type": "interrupted"}
                        )
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
                            "msg": markdown(
                                resp_text, extensions=[BlankLinkExtension()]
                            ),
                            "source": "assistant",
                        }
                        await self.send(text_data=json.dumps({"text": message}))
                        await self.send(
                            text_data=json.dumps({"action": "ongoing_turn"})
                        )
                        await self.save_message_to_db(message)

                    else:
                        all_texts += f"No text or function_call in part: {part}"
                await asyncio.sleep(0.5)

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
