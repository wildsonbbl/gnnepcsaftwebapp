"helper functions for agents"

import json
from typing import Dict, List
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlopen

from . import logger


def get_gemini_models() -> Dict[str, List[str]]:
    "get full list of gemini model"
    url = (
        "https://raw.githubusercontent.com/wildsonbbl/gnnepcsaftwebapp/"
        "refs/heads/dev/gnnmodel/data/gemini_models.json"
    )
    try:
        with urlopen(url) as ans:
            ans = ans.read().decode("utf8").rstrip()
            return json.loads(ans)
    except (HTTPError, URLError) as e:
        logger.error("Error getting gemini models from github: %s", e)
        return {}


def get_ollama_models() -> Dict[str, List[Dict[str, str]]]:
    "get full list of ollama models"
    url = "http://localhost:11434/api/tags"
    try:
        with urlopen(url) as ans:
            ans = ans.read().decode("utf8").rstrip()
            return json.loads(ans)
    except (HTTPError, URLError) as e:
        logger.error("Error getting ollama models from local server: %s", e)
        return {}


def is_api_key_valid(api_key: str) -> bool:
    "Check if the API key is valid."
    url = f"https://generativelanguage.googleapis.com/v1/models?key={quote(api_key)}"
    try:
        with urlopen(url) as ans:
            ans = ans.read().decode("utf8").rstrip()
            return True
    except (HTTPError, URLError) as e:
        logger.error("Error checking Google API key: %s", e)
        return False


def is_ollama_online() -> bool:
    "Check if the Ollama server is online."
    url = "http://localhost:11434/api/tags"
    try:
        with urlopen(url) as ans:
            ans = ans.read().decode("utf8").rstrip()
            return True
    except (HTTPError, URLError) as e:
        logger.info("Ollama is offline error: %s", e)
        return False
