#!/bin/bash

rm -rf out || exit 0;
mkdir out;

GH_REPO="@github.com/rmflight/fmcorrelationbreastca.git"

FULL_REPO="https://$GH_TOKEN$GH_REPO"

FILE=$(ls -1t *.tar.gz | head -n 1)

tar xfz $FILE

cd out
git init
git config user.name "rmflight-travis"
git config user.email "travis"
cp "../fmcorrelationbreastca/inst/doc/*.html" .

git add .
git commit -m "deployed to github pages"
git push --force --quiet $FULL_REPO master:gh-pages
