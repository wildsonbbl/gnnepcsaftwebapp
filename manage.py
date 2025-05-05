#!/usr/bin/env python
"""Django's command-line utility for administrative tasks."""
import os
import sys

# pylint: disable = W0611,C0411
import bootstrap5
import debug_toolbar
import feos
import webapp.asgi
import webapp.wsgi
import whitenoise
from decouple import config
from uvicorn import run


def main():
    """Run administrative tasks."""
    os.environ.setdefault(
        "DJANGO_SETTINGS_MODULE",
        str(config("DJANGO_SETTINGS_MODULE", default="webapp.settings")),
    )
    try:
        # pylint: disable = C0415
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc
    execute_from_command_line(sys.argv)


if __name__ == "__main__":
    from gnnmodel.models import database_compatibility

    database_compatibility()
    if "uvicorn" in sys.argv:
        run(
            "webapp.asgi:application",
            port=int(sys.argv[2]) if len(sys.argv) > 2 else 19770,
        )
    elif "daphne" in sys.argv:
        from daphne.cli import CommandLineInterface

        CommandLineInterface().run(
            [
                "webapp.asgi:application",
                "-p",
                sys.argv[2] if len(sys.argv) > 2 else "19770",
            ]
        )
    else:
        main()
