# ![Metaboverse](https://raw.githubusercontent.com/Metaboverse/Metaboverse/master/docs/content/images/png/metaboverse_banner.png)

[![codecov](https://codecov.io/gh/Metaboverse/metaboverse-cli/branch/master/graph/badge.svg)](https://codecov.io/gh/Metaboverse/metaboverse-cli)
[![Github All Releases](https://img.shields.io/github/downloads/Metaboverse/Metaboverse/total.svg)]()
[![bioRxiv preprint](https://img.shields.io/badge/bioRxiv-10.1101%2F2020.06.25.171850-BF2636)](https://www.biorxiv.org/content/10.1101/2020.06.25.171850v1)
[![DOI](https://zenodo.org/badge/269683933.svg)](https://zenodo.org/badge/latestdoi/269683933)

## Using Metaboverse
To use the Metaboverse app, please click [here](https://github.com/Metaboverse/Metaboverse)

## What does Metaboverse do?
A current draft of the manuscript describing Metaboverse can be found [here](https://github.com/Metaboverse/manuscript/blob/master/output/manuscript.pdf)

### Requirements
- Python 3.6 or greater
- An internet connection for network curation
- The Metaboverse app

### Other Notes
- This product includes color specifications and designs developed by Cynthia Brewer (http://colorbrewer.org/).

### Build `metaboverse-cli` for the first time
```
conda create -n pyinstaller
conda activate pyinstaller

conda install python=3.8
pip install pyinstaller 
pip install pandas numpy scipy scikit-learn networkx requests

pyinstaller metaboverse-cli.spec
```

### Re-use `metaboverse-cli` build environment
```
conda activate pyinstaller

conda install python=3.8
pip install pyinstaller --upgrade
pip install pandas numpy scipy scikit-learn networkx requests --upgrade

pyinstaller metaboverse-cli.spec
```

### Utilizing a custom Metaboverse network
Custom networks are archived at https://github.com/Metaboverse/Custom-Networks/curated_organisms.
Once a network is curated, input it into the Metaboverse GUI .....

#### Example for building organism network
```
metaboverse-cli-windows.exe curate --output . --organism_id M_leprae --database_source custom --force_new_curation --organism_curation_file E:\projects\Custom-Networks\curated_organisms\m_leprae\m_leprae.json
```
