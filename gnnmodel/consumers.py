"consume module"

import json

from channels.generic.websocket import AsyncWebsocketConsumer


class ChatConsumer(AsyncWebsocketConsumer):
    "Chat consumer"

    async def receive(self, text_data=None, bytes_data=None):
        assert text_data is not None
        text_data_json = json.loads(text_data)

        # The consumer ChatConsumer is synchronous while the channel layer
        # methods are asynchronous. Therefore wrap the methods in async-to-sync
        await self.channel_layer.send(  # type: ignore
            self.channel_name,
            {
                "type": "chat.message",
                "text": {"msg": text_data_json["text"], "source": "user"},
            },
        )
        # model response here
        await self.channel_layer.send(  # type: ignore
            self.channel_name,
            {
                "type": "chat.message",
                "text": {"msg": "Bot says hello", "source": "assistant"},
            },
        )

    # Handles the chat.mesage event i.e. receives messages from the channel layer
    # and sends it back to the client.
    async def chat_message(self, event):
        "send chat message"
        text = event["text"]
        await self.send(text_data=json.dumps({"text": text}))
