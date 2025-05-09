"""Consumer utils for the current chat session."""

import asyncio
import base64
import io
import json
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
from google.genai.types import Part
from markdown import markdown

from . import logger
from .agents import all_tools
from .agents_utils import get_gemini_models, get_ollama_models, is_ollama_online
from .chat_utils import (
    BlankLinkExtension,
    CustomJSONEncoder,
    docstring_to_html,
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
    "gemini-2.5-flash-preview-04-17",
    "gemini-2.5-pro-preview-05-06",
    "gemini-2.5-pro-exp-03-25",
    "gemini-2.0-flash",
    "gemini-2.0-flash-lite",
    "gemini-1.5-flash",
    "gemini-1.5-pro",
    "gemini-1.5-flash-8b",
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


async def _get_mcp_server_names_from_config() -> List[str]:
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
        logger.warning("MCP config file is not valid JSON. Cannot parse server names.")
    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.error("Error reading MCP configuration file for server names: %s", e)
    return mcp_server_names


class CurrentChatSessionConsumerUtils(CurrentChatSessionConsumer):
    "Chat Consumer Utils"

    async def load_session_data(
        self,
        session: ChatSession,
        current_tool_map: Dict[str, Any],
        valid_selected_tools: List[str],
        mcp_server_names: List[str],
    ):
        "load session data"
        self.runner, self.runner_session = await start_agent_session(
            self.session_id,
            session.model_name,
            tools=[current_tool_map[tool] for tool in valid_selected_tools],
        )

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
            "msg": markdown(
                f"**Error:** {error_text}", extensions=[BlankLinkExtension()]
            ),
            "source": "assistant",
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

    async def activate_mcp_server(self, server_name_to_activate: Optional[str]):
        "activate the servers on the config"
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
        except Exception as e:  # pylint: disable=broad-exception-caught
            error_message = (
                f"Error reading MCP configuration"
                f" file: {settings.MCP_SERVER_CONFIG}. Error: {e}"
            )
            logger.error(error_message)
            mcp_server_config = {}

        if not error_message and "mcpServers" in mcp_server_config:
            logger.debug(mcp_server_config)
            for mcpserver_name in mcp_server_config["mcpServers"]:
                if (
                    server_name_to_activate
                    and mcpserver_name != server_name_to_activate
                ):
                    continue
                mcpserver = mcp_server_config["mcpServers"][mcpserver_name]
                command = mcpserver.get("command")
                args = mcpserver.get("args")
                env = mcpserver.get("env")

                if not isinstance(command, str) or not command:
                    logger.error(
                        "Invalid or missing 'command' for MCP server: %s",
                        mcpserver_name,
                    )
                    continue
                if not isinstance(args, list):
                    logger.error(
                        "Invalid or missing 'args' for MCP server: %s. Must be a list.",
                        mcpserver_name,
                    )
                    continue

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

        return activated_tool_names, error_message

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
