#!/bin/bash

echo 'hexamers <- c('
python hexamers.py | awk 'NR > 1 {printf(",\n")} {printf("%s", $1)}'
echo ')'
