"consume module"

import json
from typing import Any, Dict

from .action_handlers import ChatConsumerHandleActions
from .consumer_utils import mcp_exit_stack


class ChatConsumer(ChatConsumerHandleActions):
    "Chat consumer"

    async def connect(self):
        """Connect to the websocket"""
        await self.accept()

        session = await self.get_or_create_last_session()
        current_tools = self.original_tools + self.mcp_tools
        current_tool_map = self.get_current_tool_map(current_tools)
        valid_selected_tools = await self.validate_and_update_tools(
            session, current_tool_map
        )
        await self.load_session_data(session, current_tool_map, valid_selected_tools)

    async def disconnect(self, code):
        await mcp_exit_stack.aclose()

    async def receive(self, text_data=None, bytes_data=None):
        assert text_data is not None
        text_data_json: Dict[str, Any] = json.loads(text_data)

        if "action" in text_data_json:
            await self.handle_actions(text_data_json)
        elif "text" in text_data_json:
            text = text_data_json["text"]
            file_info = text_data_json.get("file")

            if not text.strip() and not file_info:
                return

            await self.handle_text(text, file_info)
