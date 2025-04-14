"""utils for adk agent runner"""

from google.adk.runners import Runner
from google.adk.sessions.in_memory_session_service import InMemorySessionService

from .agents import root_agent

APP_NAME = "GNNePCSAFT Agent"
session_service = InMemorySessionService()


def start_agent_session(session_id: str):
    """Starts an agent session"""

    # Create a Session
    session = session_service.create_session(
        app_name=APP_NAME,
        user_id=session_id,
        session_id=session_id,
    )

    # Create a Runner
    runner = Runner(
        app_name=APP_NAME,
        agent=root_agent,
        session_service=session_service,
    )

    return runner, session
