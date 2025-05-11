"""Consumer utils for the current chat session."""

import asyncio
import json
import os

from google.genai.types import Content, Part
from markdown import markdown

from . import logger
from .agents_utils import is_api_key_valid
from .chat_utils import USER_ID, BlankLinkExtension
from .models import ChatSession
from .session_operations import ChatSessionsDBOperations


class ChatConsumerMessagingOperations(ChatSessionsDBOperations):
    "messaging functions"

    async def agent_to_client_messaging(self, text, file_part=None):
        """Agent to client communication"""
        try:
            if not await self.handle_api_key_validation():
                return
            parts = [Part.from_text(text=text)]
            if file_part:
                parts.append(file_part)

            await self.process_user_message(parts)

        except asyncio.CancelledError:
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

    async def handle_api_key_validation(self) -> bool:
        "handle api key validation"
        current_session: ChatSession = await self.get_or_create_session()
        model_name = current_session.model_name
        key_is_valid = True

        if model_name and model_name.lower().startswith("gemini"):
            google_api_key = os.getenv("GOOGLE_API_KEY", os.getenv("GEMINI_API_KEY"))
            key_is_valid = False
            if google_api_key:
                key_is_valid = is_api_key_valid(google_api_key)

            if not key_is_valid:
                await self.send(
                    text_data=json.dumps(
                        {
                            "action": "api_key_required",
                            "model_name": model_name,
                        }
                    )
                )
                error_message_for_log = {
                    "msg": markdown(
                        f"**System Message:** Cannot process message. "
                        f"Google API Key is missing or invalid for model `{model_name}`. "
                        "Please ensure the `GOOGLE_API_KEY` or `GEMINI_API_KEY` "
                        "environment variable is correctly set. "
                        "You can find more about Gemini API Keys "
                        "<a href='https://ai.google.dev/gemini-api/docs/api-key' "
                        "target='_blank'>here</a>.",
                        extensions=[BlankLinkExtension()],
                    ),
                    "source": "assistant",
                }
                await self.send(text_data=json.dumps({"text": error_message_for_log}))
                await self.save_message_to_db(error_message_for_log)
                await self.send(
                    text_data=json.dumps(
                        {"action": "end_turn", "type": "api_key_error"}
                    )
                )
        return key_is_valid

    async def client_to_agent_messaging(self, text, file_info_for_db=None):
        """Client to agent communication"""
        user_message = {
            "msg": markdown(text, extensions=[BlankLinkExtension()]),
            "source": "user",
        }
        if file_info_for_db:
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

    async def process_user_message(self, parts):
        "process user message"
        async for event in self.runner.run_async(
            new_message=Content(role="user", parts=parts),
            user_id=USER_ID,
            session_id=self.session_id,
            run_config=self.run_config,
        ):
            if event.interrupted:
                await self.send(
                    text_data=json.dumps({"action": "end_turn", "type": "interrupted"})
                )
                break

            all_parts = event.content and event.content.parts
            if not all_parts or event.partial:
                continue
            all_texts = ""
            for part in all_parts:
                if part.text:
                    resp_text = part.text
                    all_texts += resp_text + "\n\n"
                elif part.function_call and part.function_call.name:
                    resp_text = "**Calling:** `" + part.function_call.name + "`"
                    all_texts += resp_text + "\n\n"
                else:
                    resp_text = None
                if resp_text:
                    message = {
                        "msg": markdown(resp_text, extensions=[BlankLinkExtension()]),
                        "source": "assistant",
                    }
                    await self.send(text_data=json.dumps({"text": message}))
                    await self.send(text_data=json.dumps({"action": "ongoing_turn"}))
                    await self.save_message_to_db(message)
                else:
                    all_texts += f"No text or function_call in part: {part}"
            await asyncio.sleep(0.5)
