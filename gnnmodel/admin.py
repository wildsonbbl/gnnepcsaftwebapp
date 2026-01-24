"Django admin module."

from django.contrib import admin

from .models import ThermoMLDenData, ThermoMLVPData

admin.site.register(
    ThermoMLDenData,
    admin.ModelAdmin,
)

admin.site.register(
    ThermoMLVPData,
    admin.ModelAdmin,
)
