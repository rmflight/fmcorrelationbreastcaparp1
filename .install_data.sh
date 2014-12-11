#!/bin/bash

rm -rf fmdatabreastcaparp1 || exit 0;
git clone https://github.com/rmflight/fmdatabreastcaparp1.git

mkdir fmdatabreastcaparp1/data
wget https://dl.dropbox.com/s/89xgvje4rfl0lml/fmdatabreastcaparp1.zip
unzip fmdatabreastcaparp1.zip -d fmdatabreastcaparp1/data

Rscript -e 'library(devtools); install("fmdatabreastcaparp1")'

