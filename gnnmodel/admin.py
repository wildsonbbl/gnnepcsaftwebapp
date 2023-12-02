"Django admin module."

from django.contrib import admin

from .models import GnnepcsaftDB


class GnnepcsaftDBAdmin(admin.ModelAdmin):
    "DB table row look at admin page."
    list_display = ["inchi", "counting"]


admin.site.register(GnnepcsaftDB, GnnepcsaftDBAdmin)
