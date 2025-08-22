"consume module"

import json
from typing import Any, Dict

from . import logger
from .action_handlers import ChatConsumerHandleActions


class ChatConsumer(ChatConsumerHandleActions):
    "Chat consumer"

    async def connect(self):
        """Connect to the websocket"""
        await self.accept()

        session = await self.get_or_create_last_session()
        await self.load_session_data(
            session,
        )

    async def disconnect(self, code):
        for mcp_tool_set in self.mcp_tool_sets[::-1]:
            logger.debug("Closing MCP tool set: %s", mcp_tool_set)
            await mcp_tool_set.close()
        self.mcp_tool_sets = []

    async def receive(self, text_data=None, bytes_data=None):
        assert text_data is not None
        text_data_json: Dict[str, Any] = json.loads(text_data)

        if "action" in text_data_json:
            await self.handle_actions(text_data_json)
            logger.debug("Handled action: %s", text_data_json["action"])
        elif "text" in text_data_json:
            text = text_data_json["text"]
            file_info = text_data_json.get("file")

            if not text.strip() and not file_info:
                return

            await self.handle_text(text, file_info)
            logger.debug("Handled text input")
