#!/bin/bash/

## git remote add origin <server> ## For setting up new folder and connecting it to a remo repo

git add --all * ## Add to the index

git commit -m "$1" ## Add to the staging area

git push origin master ## send to the remote repository

echo "\nGit synced!\n"

## My username: DanJeffries
## Password: fishy zebras
