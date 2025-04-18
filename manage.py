#!/usr/bin/env python
"""Django's command-line utility for administrative tasks."""
import os
import sys

# pylint: disable = W0611,C0411
import bootstrap5
import debug_toolbar
import feos
import gunicorn
import pwa
import webapp.asgi
import webapp.wsgi
import whitenoise
from decouple import config


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
    main()
