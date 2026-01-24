"Module for django models (dbs)."

from django.db import models


class ThermoMLDenData(models.Model):
    "Table at DB for ThermoML Archive and GNNPCSAFT predicted density data."

    inchi = models.CharField(max_length=255)
    den = models.JSONField()


class ThermoMLVPData(models.Model):
    "Table at DB for ThermoML Archive and GNNPCSAFT predicted vapor pressure data."

    inchi = models.CharField(max_length=255)
    vp = models.JSONField()
