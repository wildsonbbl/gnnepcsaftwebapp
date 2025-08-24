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
    model_name = models.CharField(max_length=100, default="gemini-2.5-flash")
    selected_tools = models.JSONField(default=list)
    selected_mcp_servers = models.JSONField(default=list)

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
    """Check if the local database is compatible."""

    app_db = str(settings.BASE_DIR / "gnnepcsaft.db")
    local_db = str(settings.DB_PATH)

    try:
        logger.info("Verifying database compatibility...")

        # Se o DB local não existe, tente semear a partir do template
        if not os.path.exists(local_db):
            logger.info("Local DB not found at %s", local_db)
            if os.path.exists(app_db) and os.path.abspath(app_db) != os.path.abspath(
                local_db
            ):
                # Garante que o diretório existe
                os.makedirs(os.path.dirname(local_db) or ".", exist_ok=True)
                shutil.copyfile(app_db, local_db)
                logger.warning("Database created from template %s", app_db)
            else:
                if not os.path.exists(app_db):
                    logger.warning(
                        "Template DB not found at %s; skipping seeding.", app_db
                    )
            return

        # Não tente abrir o template se ele não existir (evita criar DB vazio)
        if not os.path.exists(app_db):
            logger.error(
                "Template DB not found at %s; skipping "
                "compatibility check to avoid creating a blank DB.",
                app_db,
            )
            return

        def fetch_example(db_path: str):
            try:
                with sqlite3.connect(db_path) as conn:
                    cur = conn.cursor()
                    row = cur.execute(
                        "SELECT den FROM gnnmodel_thermomldendata WHERE inchi = ?;",
                        ("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",),
                    ).fetchone()
                return row, None
            except sqlite3.Error as e:
                return None, e

        example_local, err_local = fetch_example(local_db)
        example_app, err_app = fetch_example(app_db)

        # Trate qualquer erro ou ausência do registro sentinela como incompatível
        incompatible = (
            err_local is not None
            or err_app is not None
            or example_local is None
            or example_app is None
            or example_local != example_app
        )

        if incompatible:
            logger.warning(
                "Incompatible database detected (local=%s, app=%s, err_local=%s, err_app=%s)",
                example_local,
                example_app,
                err_local,
                err_app,
            )
            # Substituir apenas se o template for distinto do local
            if os.path.abspath(app_db) != os.path.abspath(local_db):
                backup_path = f"{local_db}.broken_{int(time.time())}"
                try:
                    os.rename(local_db, backup_path)
                    logger.warning(
                        "Broken database detected. Backup saved in %s", backup_path
                    )
                except OSError as ex:
                    logger.error("Failed to backup DB %s: %s", local_db, ex)
                try:
                    shutil.copyfile(app_db, local_db)
                    logger.warning("Database substituted by new copy from %s", app_db)
                except OSError as ex:
                    logger.error(
                        "Failed to copy template DB from %s to %s: %s",
                        app_db,
                        local_db,
                        ex,
                    )
        else:
            logger.info("Database is compatible.")
    except Exception as ex:  # pylint: disable=broad-except
        logger.error("Error while checking DB compatibility: %s", ex)
