"consume module"

import asyncio
import json
import uuid
from datetime import datetime

from channels.db import database_sync_to_async
from channels.generic.websocket import AsyncWebsocketConsumer
from google.adk.agents.run_config import RunConfig
from google.adk.runners import Runner
from google.adk.sessions.in_memory_session_service import Session
from google.genai.types import Content, Part
from markdown import markdown
from markdown.extensions import Extension
from markdown.treeprocessors import Treeprocessor

from .chat_utils import APP_NAME, USER_ID, session_service, start_agent_session
from .models import ChatSession


class CustomJSONEncoder(json.JSONEncoder):
    """Custom JSON encoder that handles UUIDs and datetimes"""

    def default(self, o):
        if isinstance(o, uuid.UUID):
            return str(o)
        if isinstance(o, datetime):
            return o.isoformat()
        return super().default(o)


class BlankLinkExtension(Extension):
    "extension to set links to target=_blank"

    def extendMarkdown(self, md):
        md.treeprocessors.register(BlankLinkProcessor(md), "blank_link_processor", 15)


class BlankLinkProcessor(Treeprocessor):  # pylint: disable=too-few-public-methods
    "processor to set links to target=_blank"

    def run(self, root):
        for element in root.iter("a"):
            element.set("target", "_blank")
        return root


class ChatConsumer(AsyncWebsocketConsumer):
    "Chat consumer"

    session_id: str
    runner_session: Session
    runner: Runner
    run_config = RunConfig(response_modalities=["TEXT"])

    async def connect(self):
        """Connect to the websocket"""
        await self.accept()

        # Get the last session or create a new one if none exists
        last_session = await self.get_last_session()
        if last_session:
            self.session_id = str(last_session.session_id)
            session_name = last_session.name
        else:
            # Create a default session if none exists
            self.session_id = str(uuid.uuid4())
            session_name = "New Session"
            await self.create_new_session(self.session_id, session_name)

        # Initialize the agent with the session
        self.runner, self.runner_session = start_agent_session(self.session_id)

        # Send the current session info to the client
        await self.send(
            text_data=json.dumps(
                {
                    "action": "session_loaded",
                    "session_id": self.session_id,
                    "name": session_name,
                },
                cls=CustomJSONEncoder,
            )
        )

        # Load previous messages if session exists
        if last_session and last_session.messages:
            await self.send(
                text_data=json.dumps(
                    {"action": "load_messages", "messages": last_session.messages},
                    cls=CustomJSONEncoder,
                )
            )

    @database_sync_to_async
    def get_last_session(self):
        """Get the most recently updated session"""
        try:
            return ChatSession.objects.order_by("-updated_at").first()
        except ChatSession.DoesNotExist:
            return None

    @database_sync_to_async
    def create_new_session(self, session_id, name):
        """Create a new session"""
        return ChatSession.objects.create(session_id=session_id, name=name)

    @database_sync_to_async
    def get_or_create_session(self, session_id):
        """Get or create a session"""
        try:
            # Try to get existing session
            return ChatSession.objects.get(session_id=session_id)
        except ChatSession.DoesNotExist:
            # Create new session
            return ChatSession.objects.create(session_id=session_id)

    @database_sync_to_async
    def save_message_to_db(self, message):
        """Save message to database"""
        session = ChatSession.objects.get(session_id=self.session_id)
        session.add_message(message)
        return session

    def serialize_session(self, session):
        """Convert session data to JSON-serializable format"""
        return {
            "session_id": str(session["session_id"]),  # Convert UUID to string
            "name": session["name"],
            "created_at": session[
                "created_at"
            ].isoformat(),  # Convert datetime to ISO format string
            "updated_at": session["updated_at"].isoformat(),
        }

    @database_sync_to_async
    def get_all_sessions(self):
        """Get all sessions"""
        sessions = list(
            ChatSession.objects.values("session_id", "name", "created_at", "updated_at")
        )
        # Convert each session to a JSON-serializable format
        return [self.serialize_session(session) for session in sessions]

    async def receive(self, text_data=None, bytes_data=None):
        assert text_data is not None
        text_data_json = json.loads(text_data)

        # Handle different types of messages
        if "action" in text_data_json:
            action = text_data_json["action"]

            if action == "delete_session":
                session_id = text_data_json["session_id"]
                success = await self.delete_session(session_id)

                # If the deleted session was the current one, load the most recent session
                if success and session_id == self.session_id:
                    last_session = await self.get_last_session()
                    if last_session:
                        self.session_id = str(last_session.session_id)
                        session_name = last_session.name
                        self.runner, self.runner_session = start_agent_session(
                            self.session_id
                        )

                        # Send the current session info to the client
                        await self.send(
                            text_data=json.dumps(
                                {
                                    "action": "session_loaded",
                                    "session_id": self.session_id,
                                    "name": session_name,
                                },
                                cls=CustomJSONEncoder,
                            )
                        )

                        # Load previous messages
                        await self.send(
                            text_data=json.dumps(
                                {
                                    "action": "load_messages",
                                    "messages": last_session.messages,
                                },
                                cls=CustomJSONEncoder,
                            )
                        )
                    else:
                        # Create a new session if no other sessions exist
                        new_session_id = str(uuid.uuid4())
                        new_session_name = "New Session"
                        await self.create_new_session(new_session_id, new_session_name)
                        self.session_id = new_session_id
                        self.runner, self.runner_session = start_agent_session(
                            new_session_id
                        )

                        await self.send(
                            text_data=json.dumps(
                                {
                                    "action": "session_loaded",
                                    "session_id": new_session_id,
                                    "name": new_session_name,
                                },
                                cls=CustomJSONEncoder,
                            )
                        )

                        # Clear messages
                        await self.send(
                            text_data=json.dumps(
                                {"action": "load_messages", "messages": []}
                            )
                        )

                # Send confirmation and refresh sessions list
                await self.send(
                    text_data=json.dumps(
                        {
                            "action": "session_deleted",
                            "success": success,
                            "session_id": session_id,
                        }
                    )
                )

                # Send updated sessions list
                sessions = await self.get_all_sessions()
                await self.send(
                    text_data=json.dumps(
                        {"action": "sessions_list", "sessions": sessions},
                        cls=CustomJSONEncoder,
                    )
                )

            elif action == "get_sessions":
                # Return list of available sessions
                sessions = await self.get_all_sessions()

                await self.send(
                    text_data=json.dumps(
                        {"action": "sessions_list", "sessions": sessions}
                    )
                )

            elif action == "create_session":
                # Create a new session
                new_session_id = str(uuid.uuid4())
                name = text_data_json.get("name", "New Session")
                session = await database_sync_to_async(ChatSession.objects.create)(
                    session_id=new_session_id, name=name
                )

                # Reset the agent session
                self.session_id = new_session_id
                self.runner, self.runner_session = start_agent_session(new_session_id)

                await self.send(
                    text_data=json.dumps(
                        {
                            "action": "session_created",
                            "session_id": new_session_id,
                            "name": name,
                        }
                    )
                )

            elif action == "load_session":
                # Load an existing session
                session_id = text_data_json["session_id"]
                self.session_id = session_id
                self.runner, self.runner_session = start_agent_session(session_id)

                session = await self.get_or_create_session(session_id)
                await self.send(
                    text_data=json.dumps(
                        {"action": "load_messages", "messages": session.messages}
                    )
                )

            elif action == "rename_session":
                # Rename a session
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

        elif "text" in text_data_json:
            if text_data_json["text"] == "":
                return

            # Store user message in database
            user_message = {"msg": text_data_json["text"], "source": "user"}
            await self.save_message_to_db(user_message)

            await self.client_to_agent_messaging(text_data_json["text"])
            await self.agent_to_client_messaging(text_data_json["text"])

    @database_sync_to_async
    def delete_session(self, session_id):
        """Delete a session"""
        try:
            session = ChatSession.objects.get(session_id=session_id)
            session.delete()

            session_service.delete_session(
                app_name=APP_NAME, user_id=USER_ID, session_id=session_id
            )

            return True
        except ChatSession.DoesNotExist:
            return False

    async def agent_to_client_messaging(self, text):
        """Agent to client communication"""
        try:
            async for event in self.runner.run_async(
                new_message=Content(role="user", parts=[Part.from_text(text=text)]),
                user_id=USER_ID,
                session_id=self.session_id,
                run_config=self.run_config,
            ):

                if event.interrupted:
                    print("[INTERRUPTED]")
                    await self.send(text_data=json.dumps({"action": "end_turn"}))
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
                        message = {
                            "msg": markdown(text, extensions=[BlankLinkExtension()]),
                            "source": "assistant",
                        }
                        await self.send(text_data=json.dumps({"text": message}))
                        await self.send(
                            text_data=json.dumps({"action": "ongoing_turn"})
                        )
                        await self.save_message_to_db(message)

                    else:
                        all_texts += f"No text or function_call in part: {part}"
                await asyncio.sleep(0.5)
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error with agent: {e}")
            message = {
                "msg": markdown(
                    f"***Error with agent**: `{e}`*", extensions=[BlankLinkExtension()]
                ),
                "source": "assistant",
            }
            await self.send(text_data=json.dumps({"text": message}))
            await self.save_message_to_db(message)
        finally:
            await self.send(text_data=json.dumps({"action": "end_turn"}))

    async def client_to_agent_messaging(self, text):
        """Client to agent communication"""
        await self.send(
            text_data=json.dumps({"text": {"msg": text, "source": "user"}}),
        )
        await self.send(text_data=json.dumps({"action": "ongoing_turn"}))
        # print(f"[CLIENT TO AGENT]: {text}")

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
            # print(f"[AGENT TO CLIENT]: {text}")
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
        # print("Done")
