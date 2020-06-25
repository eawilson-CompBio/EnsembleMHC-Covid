#!/bin/bash

dataset=$(mhcflurry-downloads path models_class1_presentation)

mhcflurry-predict --allele-column HLA --peptide-column peptide --models "$dataset/models/" --out mhcflurry_pred.out $1

