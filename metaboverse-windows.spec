# -*- mode: python3.8 ; coding: utf-8 -*-
import sys
sys.setrecursionlimit(5000)

block_cipher = None

a = Analysis(
            ['metaboverse_cli\\__main__.py'],
             pathex=[
              'metaboverse_cli'],
             binaries=[
               ('metaboverse_cli\\analyze\\data\\metabolite_mapping.pickle.zip',
               'analyze\\data')
             ],
             datas=[
              ('README.md', '.')
             ],
             hiddenimports=[
              'scipy.special.cython_special',
              'scipy.spatial.transform._rotation_groups'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[
              'metaboverse_cli\\test\\',
              'metaboverse_cli\\__test__.py',
              'metaboverse_cli\\curate\\test\\',
              'metaboverse_cli\\curate\\__test__.py',
              'metaboverse_cli\\mapper\\test\\',
              'metaboverse_cli\\mapper\\__test__.py',
              'metaboverse_cli\\analyze\\test\\',
              'metaboverse_cli\\analyze\\__test__.py'
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
          name='metaboverse-cli-win',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )
