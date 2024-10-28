"Module for django models (dbs)."
from django.db import models


class GnnepcsaftPara(models.Model):
    "Table at DB for gnnepcsaft predicted params."
    inchi = models.CharField(max_length=255)
    m = models.FloatField()
    sigma = models.FloatField()
    e = models.FloatField()
    plot_den = models.CharField(max_length=50)
    plot_vp = models.CharField(max_length=50)
    plot_mol = models.CharField(max_length=50)
