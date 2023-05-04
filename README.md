# ![Metaboverse](https://raw.githubusercontent.com/Metaboverse/Metaboverse/master/docs/content/images/png/metaboverse_banner.png)

[![Github All Releases](https://img.shields.io/github/downloads/Metaboverse/Metaboverse/total.svg)](https://github.com/Metaboverse/Metaboverse/releases/)
[![manuscript](https://img.shields.io/badge/Manuscript-10.1038%2Fs41556--023--01117--9-red)](https://doi.org/10.1038/s41556-023-01117-9)
[![DOI](https://zenodo.org/badge/269683933.svg)](https://zenodo.org/badge/latestdoi/269683933)

## What does Metaboverse do?
Integrating multi- or single-omic metabolic data upon the metabolic network can be challenging for a variety of reasons. Metaboverse seeks to simplify this task for users by providing a simple, user-friendly interface for layering their data on a dynamic representation of the metabolic network and automatically searching the network for interesting regulatory or other patterns. Additionally, Metaboverse provides several tools to enable the contextualization of metabolic data.

## Citing Metaboverse
If you use Metaboverse in your data analysis, please cite as:
```
Berg JA, Zhou Y, Ouyang Y, Cluntun AA, Waller TC, Conway ME, Nowinski SM, Van Ry T, George I,
Cox JE, Wang B, Rutter J.
Metaboverse enables automated discovery and visualization of diverse metabolic regulatory patterns.
Nature Cell Biology. (2023) doi: https://doi.org/10.1038/s41556-023-01117-9
```

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

conda install python
pip install pyinstaller 
pip install -r requirements.txt

pyinstaller metaboverse-cli.spec
```

### Utilizing a custom Metaboverse network (in progress)
Custom networks are archived at https://github.com/Metaboverse/Custom-Networks/curated_organisms.
Once a network is curated, input it into the Metaboverse GUI .....

#### Example for building organism network
```
metaboverse-cli-windows.exe curate --output . --organism_id M_leprae --database_source custom --force_new_curation --organism_curation_file E:\projects\Custom-Networks\curated_organisms\m_leprae\m_leprae.json
```
