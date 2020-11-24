#!/bin/bash

conda create -q -n mhcflurry-env python=3.6 'keras==2.3.1'

conda activate mhcflurry-env

pip install mhcflurry==1.6.0

mhcflurry-downloads fetch models_class1_presentation
