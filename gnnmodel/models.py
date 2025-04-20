"Module for django models (dbs)."

import uuid

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


class ChatSession(models.Model):
    """Model to store chat sessions"""

    session_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=100, default="Unnamed Session")
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    messages = models.JSONField(default=list)
    model_name = models.CharField(max_length=100, default="gemini-2.0-flash")

    def __str__(self):
        return f"{self.name} ({self.session_id})"

    def add_message(self, message):
        """Add a message to the session"""
        self.messages.append(message)  # pylint: disable=E1101
        self.save()

    def get_messages(self):
        """Get all messages in the session"""
        return self.messages
