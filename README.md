# README
1. Clone the repo

2. Get expression data:
  - to get all data (adult and fetal) directly from Allen institute run:
  >./download_expression_data.sh
  - this should download zip files and unzip them into folders in the ./data/raw/ directory
  - if you download the data differently, make sure to edit the config.py file to point to where you've placed your data

3. Setup python3.6 environment with essential dependencies in requirements.txt
  - Use conda environments create an env named molecularAN:

  >conda create -n molecularAN python=3.6
  pip install -r requirements.txt

4. Process the Allen data to get the aggregated data:

  >python data_processing.py

  You should have a matrix each for the adult and fetal brains.

5. To get results tables:
  >python results_tables.py

6. Download SVGs to be used for figures:
  >python svg_download.py atlas_name

  where atlas_name is adult, fetal15, fetal21 or fetal21_brainstem

7. To generate SVGs:
  >python modify_svg.py
