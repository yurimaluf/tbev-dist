#!/bin/bash

git add .
if [ ${#1} -lt 3 ]
then
CURRENTDATE=`date +"%d-%b-%Y %H:%M"`;
git commit -m "${CURRENTDATE}"
else
git commit -m "${1}"
fi
git push origin main