"Module for django models (dbs)."

import os
import shutil
import sqlite3
import time
import uuid

from django.conf import settings
from django.db import models

from . import logger


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


class ChatSession(models.Model):
    """Model to store chat sessions"""

    session_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=100, default="Unnamed Session")
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    messages = models.JSONField(default=list)
    model_name = models.CharField(max_length=100, default="gemini-2.0-flash")
    selected_tools = models.JSONField(default=list)

    def __str__(self):
        return f"{self.name} ({self.session_id})"

    def add_message(self, message):
        """Add a message to the session"""
        self.messages.append(message)  # pylint: disable=E1101
        self.save()

    def get_messages(self):
        """Get all messages in the session"""
        return self.messages


def database_compatibility():
    """Check if the local database is compatible with the django models schema."""

    app_db = settings.BASE_DIR / "gnnepcsaft.db"
    local_db = settings.DB_PATH
    conn = sqlite3.connect(local_db)
    cursor = conn.cursor()
    # code to check the tables in mydatabase
    logger.info("Verifying database compatibility...")
    tables = cursor.execute(
        "SELECT name FROM sqlite_master WHERE type='table';"
    ).fetchall()
    conn.close()
    chat_tables = [
        "sessions",
        "app_states",
        "user_states",
        "events",
    ]

    if any((chat_table,) in tables for chat_table in chat_tables):
        # Backup do banco antigo
        backup_path = f"{local_db}.broken_{int(time.time())}"
        os.rename(local_db, backup_path)
        logger.warning("Broken database detected. Backup saved in %s", backup_path)
        # Substitui pelo novo
        shutil.copyfile(app_db, local_db)
        logger.warning("Database substituted by new copy from %s", app_db)
    else:
        logger.info("Database is compatible with django models schema.")
