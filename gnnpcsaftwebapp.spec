# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import copy_metadata

datas = [
    ("./icons.json", "."),
    ("./gnnpcsaftwebapp.db", "."),
    ("./gnnmodel/templates", "./gnnmodel/templates"),
    ("./productionfiles", "./productionfiles"),
]
datas += copy_metadata("django-bootstrap-v5")


a = Analysis(
    ["manage.py"],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=["tiktoken_ext.openai_public", "tiktoken_ext", "linkify-it-py"],
    hookspath=["./hooks"],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["polars"],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="gnnpcsaftwebapp",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="gnnpcsaftwebapp",
)
