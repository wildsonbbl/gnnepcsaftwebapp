"pywebview gui"

# pylint: disable = W0611,C0411

# gui.py
import os
import socket
import sys
import threading
import webbrowser

import bootstrap5
import feos
import webview
import whitenoise
from decouple import config
from django.core.management import execute_from_command_line
from uvicorn import run

import webapp.asgi
import webapp.wsgi

webview.settings["OPEN_EXTERNAL_LINKS_IN_BROWSER"] = False


def _find_free_port():
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("localhost", 0))
    port = s.getsockname()[1]
    s.close()
    return port


def _redirect_std_streams_to_devnull():
    """Prevent PyInstaller-hidden console apps from crashing on stdout/stderr."""
    sys.stdout = open(os.devnull, "w", encoding="utf-8")  # pylint: disable=R1732
    sys.stderr = open(os.devnull, "w", encoding="utf-8")  # pylint: disable=R1732


def start_django_server(port: int):
    """Start the ASGI app on a specific local port."""
    _redirect_std_streams_to_devnull()

    # Run the server on the dynamic port without the auto-reloader
    # execute_from_command_line(
    #     [sys.argv[0], "runserver", f"127.0.0.1:{port}", "--noreload"]
    # )

    run("webapp.asgi:application", port=port, log_config=None)


class WindowAPI:
    "window api"

    def __init__(self, port):
        self.port = port

    def open_link(self, url: str):
        "open nav link in new window"
        if self.check_url(url):
            webview.create_window(
                "GNNPCSAFT - New Window",
                url,
                width=800,
                height=600,
                resizable=True,
            )
        else:
            webbrowser.open(url)

    def check_url(self, url: str):
        "check url"
        return url.startswith(f"http://localhost:{self.port}") or url.startswith("/")


def start_app():
    """Bootstrap migrations, launch the local server, and open the desktop window."""
    port = _find_free_port()

    django_thread = threading.Thread(
        target=start_django_server,
        args=(port,),
        daemon=True,
    )
    django_thread.start()

    api = WindowAPI(port)
    webview.create_window(
        "GNNPCSAFT",
        f"http://localhost:{port}",
        width=800,
        height=600,
        maximized=True,
        js_api=api,
    )
    webview.start()


if __name__ == "__main__":
    start_app()
