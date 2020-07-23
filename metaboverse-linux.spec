# -*- mode: python3.8 ; coding: utf-8 -*-
import sys
sys.setrecursionlimit(5000)
block_cipher = None

a = Analysis(
            ['metaboverse_cli/__main__.py'],
             pathex=[
              'metaboverse_cli/'],
             binaries=[],
             datas=[],
             hiddenimports=[
              'scipy.special.cython_special'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[
              'metaboverse_cli/test/',
              'metaboverse_cli/__test__.py'
             ],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='metaboverse-cli-linux',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )
