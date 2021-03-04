#!/bin/bash


git pull origin master

LAST_RUN=$(cat last_date_run.txt)

mv todays_data/* previous_days/$LAST_RUN.png

/usr/local/bin/Rscript update_cor_data.R

cd ../

git add generate_readme_image/
git commit -m "update correlation data"

git push origin master
