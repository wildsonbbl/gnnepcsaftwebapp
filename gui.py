"pywebview gui"

# gui.py
import os
import socket
import sys
import threading

# pylint: disable = W0611,C0411
import bootstrap5
import feos
import webview
import whitenoise
from decouple import config
from django.core.management import execute_from_command_line
from uvicorn import run

import webapp.asgi
import webapp.wsgi


def _find_free_port():
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("localhost", 0))
    port = s.getsockname()[1]
    s.close()
    return port


PORT = _find_free_port()


def _start_django():
    # Set the environment variable for your settings module
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "webapp.settings")

    # In PyInstaller, standard outputs can crash hidden console apps,
    # so we redirect them if necessary or run with stdout bypass
    sys.stdout = open(os.devnull, "w", encoding="utf-8")
    sys.stderr = open(os.devnull, "w", encoding="utf-8")

    # Run the server on the dynamic port without the auto-reloader
    execute_from_command_line(
        [sys.argv[0], "runserver", f"127.0.0.1:{PORT}", "--noreload"]
    )


if __name__ == "__main__":
    # 2. Start Django in a background thread
    django_thread = threading.Thread(target=_start_django, daemon=True)
    django_thread.start()

    # 3. Launch the pywebview window pointing to localhost
    webview.create_window("GNNPCSAFT", f"http://localhost:{PORT}")
    webview.start()
