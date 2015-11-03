#!/bin/bash

rm *.html

git clone https://github.com/pachterlab/sleuth.git

rm -rf sleuth/man/*.html

ls sleuth/man/*.Rd |\
  sed 's/\.Rd//' |\
  xargs -n 1 -I % sh -c 'R CMD Rdconv --type=html -o %.html %.Rd'

mv sleuth/man/*.html .

rm -rf sleuth

ls *.html |\
  sed 's/\.html//' |\
  xargs -n 1 -I % sh -c 'echo "<p> <a href=\"doc/%.html\">%</a> </p>"' >\
  ../_includes/docs/links.html
