# README
1. Clone the repo

2. Get expression data:
  - to get all data (adult and fetal) directly from Allen institute run:
  >./download_expression_data.sh
  - this should download zip files and unzip them into folders in the ./data/raw/ directory

3. Edit the config.py file to point to the data folder (the download_expression_data.sh script will create a data subfolder in the current folder).

4. Setup python3.6 environment with essential dependencies in requirements.txt
  - Use conda environments create an env named molecularAN:
```
  conda create -n molecularAN python=3.6
  source activate molecularAN
  pip install -r requirements.txt
```
  - also install rpy2, this code might be needed for osx users (thanks to Clemens Brunner)
```
export LDFLAGS=-L/Library/Frameworks/R.framework/Resources/lib
export PATH=/usr/local/opt/llvm/bin:$PATH
conda install rpy2 --no-cache
```
5. Process the Allen data to get the aggregated data:

  >python data_processing.py

  You should have a matrix each for the adult and fetal brains.

6. To get results tables:
  >python results_tables.py

7. Download SVGs to be used for figures:
  >python svg_download.py atlas_name

  where atlas_name is adult, fetal15, fetal21 or fetal21_brainstem

8. To generate SVGs:
  >python modify_svg.py
