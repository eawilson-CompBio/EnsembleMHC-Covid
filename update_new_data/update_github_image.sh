#!/bin/bash

git pull origin master

LAST_RUN=$(cat last_date_run.txt)

mv todays_data/* previous_days/$LAST_RUN.png

Rscript update_cor_data.R

cd ../

git add update_new_data/
git commit -m "update correlation data"
git push origin master
