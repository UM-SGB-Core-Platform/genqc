#!/bin/bash

module load gcc/11.2
# This version of python includes required package 'numpy' and 'pandas' on the UofM Grex cluster.
module load python/3.10.4

pip install openpyxl

python plotQC.py