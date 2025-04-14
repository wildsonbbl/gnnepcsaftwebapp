"""
ASGI config for gnnepcsaft project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/4.2/howto/deployment/asgi/
"""

import os

from decouple import config
from django.core.asgi import get_asgi_application

os.environ.setdefault(
    "DJANGO_SETTINGS_MODULE",
    str(config("DJANGO_SETTINGS_MODULE", default="webapp.settings")),
)

# Initialize Django ASGI application early to ensure the AppRegistry
# is populated before importing code that may import ORM models.
django_asgi_app = get_asgi_application()

# pylint: disable = C0413,C0411
import gnnmodel.routing
from channels.auth import AuthMiddlewareStack
from channels.routing import ProtocolTypeRouter, URLRouter
from channels.security.websocket import AllowedHostsOriginValidator

application = ProtocolTypeRouter(
    {
        "http": django_asgi_app,
        "websocket": AllowedHostsOriginValidator(
            AuthMiddlewareStack(URLRouter(gnnmodel.routing.websocket_urlpatterns))
        ),
    }
)
