"consume module"

import asyncio
import json

from channels.generic.websocket import AsyncWebsocketConsumer
from google.genai.types import Content, Part
from markdown import markdown

from .chat_utils import start_agent_session


class ChatConsumer(AsyncWebsocketConsumer):
    "Chat consumer"

    live_events, live_request_queue = start_agent_session("session_01")

    async def receive(self, text_data=None, bytes_data=None):
        assert text_data is not None
        text_data_json = json.loads(text_data)

        agent_to_client_task = asyncio.create_task(self.agent_to_client_messaging())
        client_to_agent_task = asyncio.create_task(
            self.client_to_agent_messaging(text_data_json["text"])
        )
        await asyncio.gather(agent_to_client_task, client_to_agent_task)

    # Handles the chat.mesage event i.e. receives messages from the channel layer
    # and sends it back to the client.
    async def chat_message(self, event):
        "send chat message"
        text = event["text"]
        await self.send(text_data=json.dumps({"text": text}))

    async def agent_to_client_messaging(self):
        """Agent to client communication"""

        async for event in self.live_events:
            # turn_complete
            if event.turn_complete:
                print("[TURN COMPLETE]")
                break

            if event.interrupted:
                print("[INTERRUPTED]")
                break

            # Read the Content and its first Part
            part = event.content and event.content.parts and event.content.parts[0]
            if not part or not event.partial:
                continue

            # Get the text
            text = event.content and event.content.parts and event.content.parts[0].text
            if not text:
                continue
            assert isinstance(text, str)

            # # Send the text to the client
            await self.channel_layer.send(  # type: ignore
                self.channel_name,
                {
                    "type": "chat.message",
                    "text": {"msg": markdown(text), "source": "assistant"},
                },
            )
            print(f"[AGENT TO CLIENT]: {text}")
            await asyncio.sleep(2)

    async def client_to_agent_messaging(self, text):
        """Client to agent communication"""
        await self.channel_layer.send(  # type: ignore
            self.channel_name,
            {
                "type": "chat.message",
                "text": {"msg": text, "source": "user"},
            },
        )

        content = Content(role="user", parts=[Part.from_text(text=text)])
        self.live_request_queue.send_content(content=content)
        print(f"[CLIENT TO AGENT]: {text}")
        await asyncio.sleep(2)
