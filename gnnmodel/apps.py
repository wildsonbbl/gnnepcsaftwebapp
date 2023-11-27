"Django apps module."
from django.apps import AppConfig


class GnnmodelConfig(AppConfig):
    "Class representing a Django application and its configuration."
    default_auto_field = "django.db.models.BigAutoField"
    name = "gnnmodel"
