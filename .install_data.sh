#!/bin/bash

rm -rf fmdatabreastcaparp1 || exit 0;
git clone https://github.com/rmflight/fmdatabreastcaparp1.git

mkdir fmdatabreastcaparp1/data
wget http://downloads.figshare.com/article/public/1266451.zip
unzip 1266451.zip -d fmdatabreastcaparp1/data

Rscript -e 'library(devtools); install("fmdatabreastcaparp1")'

