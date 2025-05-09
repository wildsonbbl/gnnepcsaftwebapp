"""Consumer utils for the current chat session."""

import uuid
from typing import Optional

from channels.db import database_sync_to_async

from .agents import DEFAULT_MODEL, all_tools
from .chat_utils import APP_NAME, USER_ID, session_service
from .consumer_utils import CurrentChatSessionConsumerUtils
from .models import ChatSession


class ChatSessionsDBOperations(CurrentChatSessionConsumerUtils):
    "Session DB handling"

    async def get_or_create_last_session(self) -> ChatSession:
        "get or create the last session"
        last_session = await self.get_last_session()
        if last_session:
            self.session_id = str(last_session.session_id)
            session = last_session
        else:
            self.session_id = str(uuid.uuid4())
            session = await self.get_or_create_session(self.session_id)
        return session

    @database_sync_to_async
    def get_last_session(self) -> Optional[ChatSession]:
        """Get the most recently updated session"""
        try:
            return ChatSession.objects.order_by("-updated_at").first()
        except ChatSession.DoesNotExist:
            return None

    @database_sync_to_async
    def get_or_create_session(
        self,
        session_id,
        name="New Session",
        model_name=DEFAULT_MODEL,
        selected_tools=None,
    ) -> ChatSession:
        """Get or create a session"""
        try:
            return ChatSession.objects.get(session_id=session_id)
        except ChatSession.DoesNotExist:
            if selected_tools is None:
                selected_tools = [t.__name__ for t in all_tools]
            return ChatSession.objects.create(
                session_id=session_id,
                name=name,
                model_name=model_name,
                selected_tools=selected_tools,
            )

    @database_sync_to_async
    def save_message_to_db(self, message):
        """Save message to database"""
        session = ChatSession.objects.get(session_id=self.session_id)
        session.add_message(message)

    def serialize_session(self, session):
        """Convert session data to JSON-serializable format"""
        return {
            "session_id": str(session["session_id"]),
            "name": session["name"],
            "created_at": session["created_at"].isoformat(),
            "updated_at": session["updated_at"].isoformat(),
        }

    @database_sync_to_async
    def get_all_sessions(self):
        """Get all sessions"""
        sessions = list(
            ChatSession.objects.values("session_id", "name", "created_at", "updated_at")
        )
        return [self.serialize_session(session) for session in sessions]

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
