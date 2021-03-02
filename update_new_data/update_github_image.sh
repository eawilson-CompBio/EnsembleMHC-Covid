#!/bin/bash



mv todays_data/* previous_days/

Rscript update_cor_data.R

cd ../

git add update_new_data/
git commit -m "update correlation data"
git push origin master
