#!/bin/bash

set -e

notebook=$1
nb="${notebook%.*}"

# Excecute the notebook, converting to the latest version. I got some errors
# trying to execute and convert to ReST all in one shot, so just do the
# execution here . . .
ipython nbconvert \
    --log-level=DEBUG \
    --ExecutePreprocessor.enabled=True \
    --to=notebook \
    --output=${nb}.run \
    ${notebook} \
    && mv ${nb}.run.ipynb ${notebook}

# . . . and then convert to rst.
ipython nbconvert \
    --log-level=DEBUG \
    --to=rst \
    ${notebook}

sed -i "s/\`\`/\`/g" ${nb}.rst
