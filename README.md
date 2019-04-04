# Economo - 2018

Data pipeline for Economo et al., 2018 from Svoboda Lab.

This project presents a DataJoint pipeline design for the data accompanying the papers:
>>Michael N. Economo, Sarada Viswanathan, Bosiljka Tasic, Erhan Bas, Johan Winnubst, Vilas Menon, Lucas T. Graybuck,Thuc Nghi Nguyen, Kimberly A. Smith, Zizhen Yao, Lihua Wang, Charles R. Gerfen, Jayaram Chandrashekar, Hongkui Zeng, Loren L. Looger & Karel Svoboda. "Distinct descending motor cortex
pathways and their roles in movement" (2018) Nature (https://doi.org/10.1038/s41586-018-0642-9)

The data in original MATLAB format (.mat) have been ingested into a DataJoint data pipeline. 

The data: doi: 10.25378/janelia.7007846

## Design DataJoint data pipeline 
This repository contains the **Python 3.7** code of the DataJoint data pipeline design for this dataset, as well as scripts for data ingestions and visualization.
 
![Pipeline diagram of intracellular and extracellular](images/all_erd.png)

## Conversion to NWB 2.0
This repository contains the **Python 3.7** code to convert the DataJoint pipeline into NWB 2.0 format (See https://neurodatawithoutborders.github.io/)
Each NWB file represents one recording session. The conversion script can be found [here](scripts/datajoint_to_nwb.py)

## Demonstration of the data pipeline
Data queries and usages are demonstrated in this [Jupyter Notebook](notebooks/Economo-2018-examples.ipynb), where several figures from the paper are reproduced. 







