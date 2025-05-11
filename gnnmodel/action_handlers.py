"""Consumer utils for the current chat session."""

import asyncio
import json
import uuid
from typing import Any, Dict, List

from channels.db import database_sync_to_async
from django.conf import settings
from markdown import markdown

from . import logger
from .agents import DEFAULT_MODEL
from .agents_utils import get_ollama_models, is_ollama_online
from .chat_utils import BlankLinkExtension, CustomJSONEncoder
from .message_operations import ChatConsumerMessagingOperations
from .models import ChatSession


class ChatConsumerHandleActions(ChatConsumerMessagingOperations):
    "handle actions class"

    async def handle_actions(self, text_data_json: Dict[str, Any]):
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
            if action in [
                "get_sessions",
                "stop_generating",
                "get_mcp_config_content",
                "update_ollama_models",
            ]:
                await handler()
            else:
                await handler(text_data_json)
        else:
            logger.warning("Unknown action received: %s", action)

    async def handle_update_ollama_models(self):
        """Checks if Ollama is online and updates the list of available models."""
        self.availabel_models = self.gemini_models.copy()
        if is_ollama_online():
            ollama_models_data = get_ollama_models()
            if (
                ollama_models_data
                and "models" in ollama_models_data
                and ollama_models_data.get("models")
            ):
                ollama_model_names = [
                    f"ollama_chat/{model['name']}"
                    for model in ollama_models_data["models"]
                ]
                self.availabel_models.extend(ollama_model_names)
        else:
            await self.send(text_data=json.dumps({"action": "ollama_offline"}))

        await self.send(
            text_data=json.dumps(
                {
                    "action": "available_models_updated",
                    "available_models": self.availabel_models,
                }
            )
        )

    async def handle_activate_mcp(self, text_data_json: Dict[str, Any]):
        "handle activate mcp action"

        server_names_list: List[str] = text_data_json.get(
            "server_names_list", ["no server name found"]
        )

        activated_tool_names, error_message = await self.activate_mcp_server(
            servers_to_process=server_names_list
        )

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
            session: ChatSession = await self.get_or_create_session()
            session.selected_mcp_servers = server_names_list
            await database_sync_to_async(session.save)(
                update_fields=["selected_mcp_servers"]
            )

            current_tool_map, valid_selected_tools = await self.start_agent_session(
                session
            )
            await self.send(
                text_data=json.dumps(
                    {
                        "action": "mcp_activated",
                        "activated_tools": activated_tool_names,
                        "available_tools": list(current_tool_map),
                        "selected_tools": valid_selected_tools,
                        "tool_descriptions": self.tool_descriptions,
                        "selected_mcp_servers": session.selected_mcp_servers,
                    }
                )
            )

    async def handle_change_tools(self, text_data_json):
        "handle change tools"
        await database_sync_to_async(
            ChatSession.objects.filter(session_id=self.session_id).update
        )(selected_tools=text_data_json["tools"])
        session: ChatSession = await self.get_or_create_session()
        current_tool_map, valid_selected_tools = await self.start_agent_session(session)
        await self.send(
            text_data=json.dumps(
                {
                    "action": "tools_changed",
                    "selected_tools": valid_selected_tools,
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
        session: ChatSession = await self.get_or_create_session()
        await self.load_session_data(
            session,
        )

    async def handle_create_session(self, text_data_json: Dict[str, Any]):
        "handle create session"
        self.session_id = str(uuid.uuid4())
        name = text_data_json.get("name", "New Session")
        model_name = text_data_json.get("model_name", DEFAULT_MODEL)
        selected_tools_names = text_data_json.get("tools", [])
        session = await self.get_or_create_session(
            name=name,
            model_name=model_name,
            selected_tools=selected_tools_names,
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
        await self.load_session_data(
            session,
        )

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
        if success and session_id == self.session_id:
            session = await self.get_or_create_last_session()
            await self.load_session_data(
                session,
            )
        await self.send(
            text_data=json.dumps(
                {
                    "action": "session_deleted",
                    "success": success,
                    "session_id": session_id,
                }
            )
        )
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
        if model_name in self.availabel_models:
            await database_sync_to_async(
                ChatSession.objects.filter(session_id=self.session_id).update
            )(model_name=model_name)
            session: ChatSession = await self.get_or_create_session()
            await self.start_agent_session(session)
            await self.send(
                text_data=json.dumps(
                    {
                        "action": "model_changed",
                        "model_name": model_name,
                    }
                )
            )
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
        processed_file_info_for_db = None
        if file_info:
            try:
                file_part, processed_file_info_for_db = (
                    await self.process_uploaded_file(file_info)
                )
            except Exception as e:  # pylint: disable=broad-exception-caught
                logger.error("Error processing file: %s", e)
                await self.send_error_message(
                    f"Failed to process file: {file_info.get('name', 'Unknown')}. Error: {e}"
                )
                return

        await self.client_to_agent_messaging(text, processed_file_info_for_db)
        if self.agent_task and not self.agent_task.done():
            self.agent_task.cancel()
        self.agent_task = asyncio.create_task(
            self.agent_to_client_messaging(text, file_part)
        )

    async def handle_get_mcp_config_content(self):
        """Handles request to get MCP server configuration content."""
        content = ""
        error_message = None
        mcp_server_names = await self._get_mcp_server_names_from_config()
        try:
            with open(settings.MCP_SERVER_CONFIG, "r", encoding="utf-8") as f:
                content = f.read()
            if content and not mcp_server_names:
                try:
                    json.loads(content)
                except json.JSONDecodeError:
                    logger.warning(
                        "MCP config file is not valid JSON. Cannot parse server names."
                    )
                    error_message = (
                        "Configuration file is not valid JSON. "
                        "Server names could not be parsed."
                    )

        except FileNotFoundError:
            error_message = "MCP configuration file not found."
            logger.warning(error_message)

        except Exception as e:  # pylint: disable=broad-exception-caught
            error_message = f"Error reading MCP configuration file: {e}"
            logger.error(error_message)

        await self.send(
            text_data=json.dumps(
                {
                    "action": "mcp_config_content",
                    "content": content,
                    "error": error_message,
                    "mcp_server_names": mcp_server_names,
                }
            )
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
            json.loads(new_content)
            with open(settings.MCP_SERVER_CONFIG, "w", encoding="utf-8") as f:
                f.write(new_content)
            logger.info(
                "MCP configuration file updated: %s", settings.MCP_SERVER_CONFIG
            )
            mcp_server_names_after_save = await self._get_mcp_server_names_from_config()
            await self.send(
                text_data=json.dumps(
                    {
                        "action": "mcp_config_saved",
                        "success": True,
                        "mcp_server_names": mcp_server_names_after_save,
                    }
                )
            )
        except json.JSONDecodeError as e:
            logger.error("Error decoding new MCP config content as JSON: %s", e)
            await self.send(
                text_data=json.dumps(
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
                text_data=json.dumps(
                    {
                        "action": "mcp_config_saved",
                        "success": False,
                        "error": str(e),
                    }
                )
            )
