"Module for django models (dbs)."
from django.db import models


class GnnepcsaftDB(models.Model):
    "Table at DB for gnnepcsaft predicted params."
    inchi = models.CharField(max_length=255)
    m = models.FloatField()
    sigma = models.FloatField()
    e = models.FloatField()
    counting = models.IntegerField()
