"consume module"

import asyncio
import json

from channels.generic.websocket import AsyncWebsocketConsumer
from channels.layers import InMemoryChannelLayer
from google.adk.agents.run_config import RunConfig
from google.genai.types import Content, FunctionCall, Part
from markdown import markdown

from .chat_utils import start_agent_session


class ChatConsumer(AsyncWebsocketConsumer):
    "Chat consumer"

    session_id = "session_01"
    runner, runner_session = start_agent_session(session_id)
    channel_layer: InMemoryChannelLayer

    # Set response modality = TEXT
    run_config = RunConfig(response_modalities=["TEXT"])

    async def receive(self, text_data=None, bytes_data=None):
        assert text_data is not None
        text_data_json = json.loads(text_data)

        await self.client_to_agent_messaging(text_data_json["text"])
        await self.agent_to_client_messaging(text_data_json["text"])

    # Handles the chat.mesage event i.e. receives messages from the channel layer
    # and sends it back to the client.
    async def chat_message(self, event):
        "send chat message"
        text = event["text"]
        await self.send(text_data=json.dumps({"text": text}))

    async def agent_to_client_messaging(self, text):
        """Agent to client communication"""

        async for event in self.runner.run_async(
            new_message=Content(role="user", parts=[Part.from_text(text=text)]),
            user_id=self.session_id,
            session_id=self.session_id,
            run_config=self.run_config,
        ):
            # turn_complete
            if event.turn_complete:
                print("[TURN COMPLETE]")
                break

            if event.interrupted:
                print("[INTERRUPTED]")
                break

            # Read the Content and its first Part
            part = event.content and event.content.parts and event.content.parts[0]
            if not part or event.partial:
                continue

            # Get the text
            text = event.content and event.content.parts and event.content.parts[0].text
            if not text:
                text = (
                    event.content
                    and event.content.parts
                    and event.content.parts[0].function_call
                )
                if isinstance(text, FunctionCall) and text.name:
                    text = "**Calling: " + text.name + " **"
                else:
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
        print(f"[CLIENT TO AGENT]: {text}")
