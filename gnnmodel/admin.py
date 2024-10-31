"Django admin module."

from django.contrib import admin

from .models import GnnepcsaftPara


class GnnepcsaftParaAdmin(admin.ModelAdmin):
    "DB table row look at admin page."
    list_display = ["m", "inchi"]


admin.site.register(GnnepcsaftPara, GnnepcsaftParaAdmin)
