#!/bin/bash

set -e

notebook=$1
jupyter nbconvert --log-level=DEBUG --execute --to=html ${notebook}
