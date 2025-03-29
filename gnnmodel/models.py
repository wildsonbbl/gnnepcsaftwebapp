"Module for django models (dbs)."

from django.db import models


class GnnepcsaftPara(models.Model):
    "Table at DB for gnnepcsaft predicted params."

    inchi = models.CharField(max_length=255)
    smiles = models.CharField(max_length=255)
    m = models.FloatField()
    sigma = models.FloatField()
    e = models.FloatField()
    k_ab = models.FloatField()
    e_ab = models.FloatField()
    mu = models.FloatField()
    na = models.IntegerField()
    nb = models.IntegerField()


class ThermoMLDenData(models.Model):
    "Table at DB for ThermoML Archive and GNNePCSAFT predicted density data."

    inchi = models.CharField(max_length=255)
    den = models.JSONField()


class ThermoMLVPData(models.Model):
    "Table at DB for ThermoML Archive and GNNePCSAFT predicted vapor pressure data."

    inchi = models.CharField(max_length=255)
    vp = models.JSONField()


def db_update(pred, inchi):
    "Updates the gnnepcsaft db."

    new_comp = GnnepcsaftPara(
        inchi=inchi,
        m=pred[0],
        sigma=pred[1],
        e=pred[2],
        k_ab=pred[3],
        e_ab=pred[4],
        mu=pred[5],
        na=pred[6],
        nb=pred[7],
    )
    new_comp.save()
    return [new_comp]
