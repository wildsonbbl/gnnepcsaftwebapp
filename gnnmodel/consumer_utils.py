"""Consumer utils for the current chat session."""

import asyncio
import base64
import io
import json
from atexit import register
from contextlib import AsyncExitStack
from typing import Any, Dict, List, Optional, Tuple

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
    SseServerParams,
    StdioServerParameters,
)
from google.genai.types import Part

from . import logger
from .agents import all_tools
from .agents_utils import get_gemini_models, get_ollama_models, is_ollama_online
from .chat_utils import (
    CustomJSONEncoder,
    docstring_to_html,
    markdown_to_html,
    start_agent_session,
)
from .models import ChatSession

# MCPToolset exit stack
mcp_exit_stack = AsyncExitStack()


def _cleanup_mcp_resources_on_exit():
    """Função síncrona para limpar recursos de mcp_exit_stack na saída do programa."""
    try:
        asyncio.run(mcp_exit_stack.aclose())
        logger.info("mcp_exit_stack fechado com sucesso via atexit.")
    except RuntimeError as e:
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

available_models = []

GEMINI_MODELS = [
    "gemini-2.5-flash",
    "gemini-2.5-pro",
    "gemini-2.0-flash",
]

gemini_models_data = get_gemini_models()
if (
    gemini_models_data
    and "models" in gemini_models_data
    and gemini_models_data["models"]
):
    GEMINI_MODELS = gemini_models_data["models"][:]

available_models.extend(GEMINI_MODELS)

if is_ollama_online():
    ollama_models_data = get_ollama_models()
    if ollama_models_data and "models" in ollama_models_data:
        ollama_model_names = [
            f"ollama_chat/{model['name']}" for model in ollama_models_data["models"]
        ]
        available_models.extend(ollama_model_names)


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
    availabel_models = available_models.copy()
    gemini_models = GEMINI_MODELS.copy()
    mcp_exit_stack = mcp_exit_stack

    async def _get_mcp_server_names_from_config(self) -> List[str]:
        """Reads MCP server names from the configuration file."""
        mcp_server_names = []
        try:
            with open(settings.MCP_SERVER_CONFIG, "r", encoding="utf-8") as f:
                content = f.read()
            mcp_config: Dict[str, Dict[str, Any]] = json.loads(content)
            if isinstance(mcp_config.get("mcpServers"), dict):
                mcp_server_names = list(mcp_config["mcpServers"].keys())
        except FileNotFoundError:
            logger.warning(
                "MCP configuration file not found when trying to read server names."
            )
        except json.JSONDecodeError:
            logger.warning(
                "MCP config file is not valid JSON. Cannot parse server names."
            )
        except Exception as e:  # pylint: disable=broad-exception-caught
            logger.error("Error reading MCP configuration file for server names: %s", e)
        return mcp_server_names


class CurrentChatSessionConsumerUtils(CurrentChatSessionConsumer):
    "Chat Consumer Utils"

    async def load_session_data(
        self,
        session: ChatSession,
    ):
        "load session data"
        await self.activate_mcp_server(servers_to_process=session.selected_mcp_servers)
        current_tool_map, valid_selected_tools = await self.start_agent_session(session)
        mcp_server_names = await self._get_mcp_server_names_from_config()

        await self.send(
            text_data=json.dumps(
                {
                    "action": "session_loaded",
                    "session_id": self.session_id,
                    "name": session.name,
                    "model_name": session.model_name,
                    "available_models": self.availabel_models,
                    "available_tools": list(current_tool_map),
                    "selected_tools": valid_selected_tools,
                    "tool_descriptions": self.tool_descriptions,
                    "mcp_config_path": str(settings.MCP_SERVER_CONFIG),
                    "mcp_server_names": mcp_server_names,
                    "selected_mcp_servers": session.selected_mcp_servers,
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

    async def start_agent_session(self, session: ChatSession):
        "Starts agent session with updated tools"
        current_tools = self.original_tools + self.mcp_tools
        current_tool_map = self.get_current_tool_map(current_tools)
        valid_selected_tools = await self.validate_and_update_tools(
            session, current_tool_map
        )
        self.runner, self.runner_session = await start_agent_session(
            self.session_id,
            session.model_name,
            tools=[current_tool_map[tool] for tool in valid_selected_tools],
        )

        return current_tool_map, valid_selected_tools

    async def validate_and_update_tools(
        self, session: ChatSession, current_tool_map: Dict[str, Any]
    ) -> List[str]:
        "validate and update tools"
        valid_selected_tools = [
            tool for tool in session.selected_tools if tool in current_tool_map
        ]
        if len(valid_selected_tools) != len(session.selected_tools):
            await database_sync_to_async(
                ChatSession.objects.filter(session_id=self.session_id).update
            )(selected_tools=valid_selected_tools)
        return valid_selected_tools

    def get_current_tool_map(self, current_tools: List[Any]) -> Dict[str, Any]:
        """get current tool map"""
        current_tool_map = {}
        for t in current_tools:
            if hasattr(t, "__name__"):
                current_tool_map[t.__name__] = t
            elif hasattr(t, "name"):
                current_tool_map[t.name] = t
        return current_tool_map

    async def send_error_message(self, error_text):
        """Sends an error message formatted as a system message."""
        error_message = {
            "msg": markdown_to_html(f"**Error:** {error_text}"),
            "source": "assistant",
        }
        await self.send(text_data=json.dumps({"text": error_message}))

    async def get_mcp(
        self,
        parameters: Tuple[
            Optional[str],
            Optional[List[str]],
            Optional[Dict[str, str]],
            Optional[str],
            Optional[Dict[str, Any]],
        ],
    ):
        "get mcp toolset"
        command, args, env, url, headers = parameters
        if command and args:
            connection_params = StdioServerParameters(
                command=command, args=args, env=env
            )
        elif url:
            connection_params = SseServerParams(url=url, headers=headers)
        else:
            return []
        tools, _ = await MCPToolset.from_server(
            connection_params=connection_params,
            async_exit_stack=self.mcp_exit_stack,
        )
        return tools

    async def activate_mcp_server(self, servers_to_process: List[str]):
        "activate the servers on the config"
        self.mcp_tools = []
        activated_tool_names = []
        error_message = None
        await self.mcp_exit_stack.aclose()
        mcp_server_config: Dict[str, Dict[str, Dict[str, Any]]] = {}

        try:
            with open(settings.MCP_SERVER_CONFIG, "r", encoding="utf-8") as config_file:
                mcp_server_config = json.load(config_file)
        except FileNotFoundError:
            error_message = "MCP configuration file not found, config file first."
            logger.error(error_message)
        except json.JSONDecodeError as e:
            error_message = (
                f"Error decoding JSON from MCP configuration"
                f" file: {settings.MCP_SERVER_CONFIG}. Error: {e}"
            )
            logger.error(error_message)
        except Exception as e:  # pylint: disable=broad-exception-caught
            error_message = (
                f"Error reading MCP configuration"
                f" file: {settings.MCP_SERVER_CONFIG}. Error: {e}"
            )
            logger.error(error_message)

        if not error_message and "mcpServers" in mcp_server_config:
            logger.debug(mcp_server_config)
            for mcpserver_name in mcp_server_config["mcpServers"]:
                if mcpserver_name not in servers_to_process:
                    continue
                mcpserver = mcp_server_config["mcpServers"][mcpserver_name]
                parameters = (
                    mcpserver.get("command"),
                    mcpserver.get("args"),
                    mcpserver.get("env"),
                    mcpserver.get("url"),
                    mcpserver.get("headers"),
                )

                parameters_valid = self.validate_mcp_parameters(
                    parameters=(mcpserver_name, *parameters)
                )
                if not parameters_valid:
                    continue

                try:
                    new_tools = await self.get_mcp(parameters=parameters)
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

        return activated_tool_names, error_message

    def validate_mcp_parameters(self, parameters):
        "validate mcp parameters"
        mcpserver_name, command, args, env, url, headers = parameters
        parameters_valid = True
        if command and not isinstance(command, str):
            logger.error(
                "Invalid or missing 'command' for MCP server: %s",
                mcpserver_name,
            )
            parameters_valid = False
        if args and not isinstance(args, List):
            logger.error(
                "Invalid or missing 'args' for MCP server: %s. Must be a list.",
                mcpserver_name,
            )
            parameters_valid = False
        if env and not isinstance(env, Dict):
            logger.error(
                "Invalid or missing 'env' for MCP server: %s. Must be a dict.",
                mcpserver_name,
            )
            parameters_valid = False
        if url and not isinstance(url, str):
            logger.error("Invalid or missing 'url' for MCP server: %s", mcpserver_name)
            parameters_valid = False
        if headers and not isinstance(headers, Dict):
            logger.error(
                "Invalid or missing 'headers' for MCP server: %s. Must be a dict.",
                mcpserver_name,
            )
            parameters_valid = False
        return parameters_valid

    async def process_uploaded_file(self, file_info: dict):
        """Processes base64 encoded file data into a Gemini Part."""
        try:
            header, encoded_data = file_info["data"].split(",", 1)
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
                "data": file_info["data"] if mime_type.startswith("image/") else None,
            }

            if mime_type == "application/pdf":
                pdf_reader = PyPDF2.PdfReader(io.BytesIO(file_bytes))
                pdf_text = f"Content from PDF '{file_name}':\n"
                for page in pdf_reader.pages:
                    pdf_text += page.extract_text() + "\n"
                return Part.from_text(text=pdf_text), processed_file_info_for_db
            if mime_type.startswith("image/"):
                return (
                    Part.from_bytes(data=file_bytes, mime_type=mime_type),
                    processed_file_info_for_db,
                )
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
            raise
