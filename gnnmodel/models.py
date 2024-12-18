"Module for django models (dbs)."
from django.db import models


class GnnepcsaftPara(models.Model):
    "Table at DB for gnnepcsaft predicted params."
    inchi = models.CharField(max_length=255)
    smiles = models.CharField(max_length=255)
    m = models.FloatField()
    sigma = models.FloatField()
    e = models.FloatField()


class ThermoMLDenData(models.Model):
    "Table at DB for ThermoML Archive and GNNePCSAFT predicted density data."
    inchi = models.CharField(max_length=255)
    T = models.FloatField()  # Temperature
    den_tml = models.FloatField()
    den_gnn = models.FloatField(null=True)
    den_ra = models.FloatField(null=True)


class ThermoMLVPData(models.Model):
    "Table at DB for ThermoML Archive and GNNePCSAFT predicted vapor pressure data."
    inchi = models.CharField(max_length=255)
    T = models.FloatField()  # Temperature
    vp_tml = models.FloatField()
    vp_gnn = models.FloatField(null=True)
    vp_ra = models.FloatField(null=True)


def db_update(pred, inchi):
    "Updates the gnnepcsaft db."

    new_comp = GnnepcsaftPara(
        inchi=inchi,
        m=pred[0],
        sigma=pred[1],
        e=pred[2],
    )
    new_comp.save()
    return [new_comp]
