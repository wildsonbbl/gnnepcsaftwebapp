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


def db_update(pred, inchi, comp, plots):
    "Updates the gnnepcsaft db."
    if len(comp) == 0:
        new_comp = GnnepcsaftPara(
            inchi=inchi,
            m=pred[0],
            sigma=pred[1],
            e=pred[2],
            plot_den=plots[0],
            plot_vp=plots[1],
            plot_mol=plots[2],
        )
        new_comp.save()
        return [new_comp]
