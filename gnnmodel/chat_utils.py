"""utils for adk agent runner"""

from django.conf import settings
from google.adk.runners import Runner
from google.adk.sessions.database_session_service import DatabaseSessionService

from .agents import root_agent

APP_NAME = "GNNePCSAFT Agent"
db_path = settings.BASE_DIR / "gnnepcsaft.db"
DB_URL = "sqlite:///" + str(db_path)
session_service = DatabaseSessionService(DB_URL)
USER_ID = "LOCAL_USER_01"


def get_sessions_ids():
    """Get the list of sessions ids"""
    active_sessions = session_service.list_sessions(
        app_name=APP_NAME, user_id=USER_ID
    ).model_dump()
    sessions_ids = [session["id"] for session in active_sessions["sessions"]]
    return sessions_ids


def start_agent_session(session_id: str):
    """Starts an agent session"""
    sessions_ids = get_sessions_ids()

    # Create a Session
    if session_id not in sessions_ids:
        session = session_service.create_session(
            app_name=APP_NAME,
            user_id=USER_ID,
            session_id=session_id,
        )
    else:
        session = session_service.get_session(
            app_name=APP_NAME,
            user_id=USER_ID,
            session_id=session_id,
        )
        assert session

    # Create a Runner
    runner = Runner(
        app_name=APP_NAME,
        agent=root_agent,
        session_service=session_service,
    )

    return runner, session
