"Django admin module."

from django.contrib import admin

from .models import GnnepcsaftPara, ThermoMLDenData, ThermoMLVPData


class TablesAdmin(admin.ModelAdmin):
    "DB table row look at admin page."
    list_display = ["smiles", "inchi"]


admin.site.register(
    GnnepcsaftPara,
    TablesAdmin,
)

admin.site.register(
    ThermoMLDenData,
    admin.ModelAdmin,
)

admin.site.register(
    ThermoMLVPData,
    admin.ModelAdmin,
)
