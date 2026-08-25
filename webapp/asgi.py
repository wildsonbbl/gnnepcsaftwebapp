"""
ASGI config for gnnpcsaftwebapp.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/4.2/howto/deployment/asgi/
"""

import os
import sys

from decouple import config
from django.core.asgi import get_asgi_application

from .cpu_compat import show_compatibility_warning, supports_avx2

if not os.environ.get("GNNPCSAFTWEBAPP_RTCOMPAT") and not supports_avx2():
    show_compatibility_warning()
    sys.exit(1)


os.environ.setdefault(
    "DJANGO_SETTINGS_MODULE",
    str(config("DJANGO_SETTINGS_MODULE", default="webapp.settings")),
)

# Initialize Django ASGI application early to ensure the AppRegistry
# is populated before importing code that may import ORM models.
application = get_asgi_application()
