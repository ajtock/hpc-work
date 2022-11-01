#!/bin/bash

# Created: 02/08/2021
# By: Andy Tock (ajt200@cam.ac.uk)

# Create a self-contained conda environment for
# installation of python_3.9.6 and R 

mamba env create --name python_3.9.6 \
                 --file python_3.9.6_environment.yaml
