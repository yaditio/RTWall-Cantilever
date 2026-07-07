# -*- mode: python ; coding: utf-8 -*-
import os
import sys
import glob
from PyInstaller.utils.hooks import collect_submodules

block_cipher = None

# Dynamically locate the current Conda environment prefix
conda_prefix = sys.prefix
library_bin = os.path.join(conda_prefix, "Library", "bin")

# Find and include all Conda DLLs (like cairo, glib, libpng, zlib, openblas, etc.)
# so that the app is self-contained and does not require Conda or Cairo on other machines.
dll_files = glob.glob(os.path.join(library_bin, "*.dll"))
binaries = []
for dll in dll_files:
    binaries.append((dll, '.'))

# Programmatically collect all submodules for the heavier packages to ensure they are captured
hidden_imports_list = [
    'openseespy',
    'openseespywin',
    'openseespywin.opensees',
    'streamlit',
    'matplotlib',
    'numpy',
    'pandas',
]
for pkg in ['openpile', 'groundhog', 'sectionproperties', 'concreteproperties', 'opsvis', 'drawsvg']:
    hidden_imports_list += collect_submodules(pkg)

a = Analysis(
    ['run_app.py'],
    pathex=[],
    binaries=binaries,
    datas=[
        ('App.py', '.'),
    ],
    hiddenimports=hidden_imports_list,
    hookspath=['./hooks'],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='RTWallCantilever',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,  # Set to True so console output is visible for troubleshooting
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
