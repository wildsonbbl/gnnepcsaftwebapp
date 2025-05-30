"""utils for adk agent runner"""

import json
import re
import textwrap
import uuid
from datetime import datetime
from typing import Callable, List, Optional

from django.conf import settings
from google.adk.artifacts.in_memory_artifact_service import InMemoryArtifactService
from google.adk.runners import Runner
from google.adk.sessions.database_session_service import DatabaseSessionService
from markdown_it import MarkdownIt
from mdit_py_plugins.dollarmath import dollarmath_plugin

from .agents import DEFAULT_MODEL, create_root_agent

md_parser = MarkdownIt("gfm-like")
md_parser.use(dollarmath_plugin)

md_parser.render(r"$$J_w = \frac{Q}{ADT}$$")

APP_NAME = "GNNePCSAFT_Agent"
DB_URL = "sqlite:///" + str(settings.DB_CHAT_PATH)
session_service = DatabaseSessionService(DB_URL)
USER_ID = "LOCAL_USER_01"
artifact_service = InMemoryArtifactService()


def get_sessions_ids():
    """Get the list of sessions ids"""
    active_sessions = session_service.list_sessions(
        app_name=APP_NAME, user_id=USER_ID
    ).model_dump()
    sessions_ids = [session["id"] for session in active_sessions["sessions"]]
    return sessions_ids


async def start_agent_session(
    session_id: str,
    model_name: str = DEFAULT_MODEL,
    tools: Optional[List[Callable]] = None,
):
    """Starts an agent session"""
    sessions_ids = get_sessions_ids()

    root_agent = await create_root_agent(model_name, tools)

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
        artifact_service=artifact_service,
    )

    return runner, session


def docstring_to_html(doc):
    """Convert docstring to HTML"""
    if not doc:
        return ""
    # Remove indentação comum
    doc = textwrap.dedent(doc).strip()
    # Remove linhas com 4+ espaços no início (exceto dentro de crases triplas)
    lines = doc.splitlines()
    in_code_block = False
    cleaned_lines = []
    for line in lines:
        if line.strip().startswith("```"):
            in_code_block = not in_code_block
        if not in_code_block:
            # Remove indentação de 4+ espaços
            line = re.sub(r"^ {4,}", "", line)
        cleaned_lines.append(line)
    doc = "\n".join(cleaned_lines)
    # Opcional: transformar "Args:" em <b>Args:</b> para melhor visualização
    doc = re.sub(r"^Args:", "<b>Args:</b>", doc, flags=re.MULTILINE)
    # Opcional: transformar listas de argumentos em listas HTML
    doc = re.sub(
        r"^\s{0,4}(\w+)\s*\(([^)]+)\):",
        r"<li><b>\1</b> (<code>\2</code>):",
        doc,
        flags=re.MULTILINE,
    )
    doc = re.sub(
        r"^\s{0,4}(\w+)\s*:",
        r"<li><b>\1</b>:",
        doc,
        flags=re.MULTILINE,
    )
    doc = re.sub(r"`([^`]+)`", r"<code>\1</code>", doc)
    # Fechar listas HTML se houver
    if "<li>" in doc:
        doc = doc.replace("<b>Args:</b>", "<b>Args:</b><ul>")
        doc += "</ul>"
    # Agora sim, passar pelo markdown para processar crases e parágrafos
    return md_parser.render(doc)


class CustomJSONEncoder(json.JSONEncoder):
    """Custom JSON encoder that handles UUIDs and datetimes"""

    def default(self, o):
        if isinstance(o, uuid.UUID):
            return str(o)
        if isinstance(o, datetime):
            return o.isoformat()
        return super().default(o)


def markdown_to_html(txt: str):
    """Convert markdown to HTML"""
    if not txt:
        return ""
    txt = textwrap.dedent(txt).strip()
    txt = txt.replace("<think>", "<think>\n\n").replace("</think>", "\n\n</think>")

    return md_parser.render(txt)
