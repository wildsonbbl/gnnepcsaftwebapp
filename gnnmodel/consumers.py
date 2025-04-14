"consume module"

import asyncio
import json

from channels.generic.websocket import AsyncWebsocketConsumer
from channels.layers import InMemoryChannelLayer
from google.adk.agents.run_config import RunConfig
from google.genai.types import Content, Part
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
        if text_data_json["text"] == "":
            return

        await self.client_to_agent_messaging(text_data_json["text"])
        await self.agent_to_client_messaging(text_data_json["text"])
        # await self.bot_to_client_messaging()

    async def agent_to_client_messaging(self, text):
        """Agent to client communication"""

        async for event in self.runner.run_async(
            new_message=Content(role="user", parts=[Part.from_text(text=text)]),
            user_id=self.session_id,
            session_id=self.session_id,
            run_config=self.run_config,
        ):

            if event.interrupted:
                print("[INTERRUPTED]")
                await self.send(
                    text_data=json.dumps(
                        {
                            "text": {
                                "msg": "Turn interrupted, brother",
                                "source": "turn_end",
                                "end_turn": True,
                            }
                        }
                    ),
                )
                break

            all_parts = event.content and event.content.parts
            if not all_parts or event.partial:
                continue
            all_texts = ""
            for part in all_parts:
                if part.text:
                    text = part.text
                    all_texts += text + "\n\n"
                elif part.function_call and part.function_call.name:
                    text = "**Calling:** `" + part.function_call.name + "`"
                    all_texts += text + "\n\n"
                else:
                    text = None
                if text:
                    await self.send(
                        text_data=json.dumps(
                            {
                                "text": {
                                    "msg": markdown(text),
                                    "source": "assistant",
                                    "end_turn": False,
                                }
                            }
                        ),
                    )
                else:
                    all_texts += f"No text or function_call in part: {part}"
            print(f"[AGENT TO CLIENT]: {all_texts}".encode("utf-8"))
            await asyncio.sleep(0.5)
        print("[TURN COMPLETE]")
        await self.send(
            text_data=json.dumps(
                {
                    "text": {
                        "msg": "Turn completed, brother",
                        "source": "turn_end",
                        "end_turn": True,
                    }
                }
            ),
        )

    async def client_to_agent_messaging(self, text):
        """Client to agent communication"""
        await self.send(
            text_data=json.dumps(
                {"text": {"msg": text, "source": "user", "end_turn": False}}
            ),
        )
        print(f"[CLIENT TO AGENT]: {text}")

    async def bot_to_client_messaging(self):
        """Bot to client communication, for testing"""
        textes = ["booot mf", "booot mf 2", "booot mf 3", "'I'M A BOT MF"]
        for text in textes:
            await self.send(
                text_data=json.dumps(
                    {"text": {"msg": text, "source": "assistant", "end_turn": False}}
                ),
            )
            print(f"[AGENT TO CLIENT]: {text}")
            await asyncio.sleep(2)
        await self.send(
            text_data=json.dumps(
                {
                    "text": {
                        "msg": "Turn completed, brother",
                        "source": "turn_end",
                        "end_turn": True,
                    }
                }
            ),
        )
        print("Done")
