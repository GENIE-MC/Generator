#!/usr/bin/bash

shopt -s nullglob
for f in *.*
do
    echo "file - $f"
    sed -f `echo $GENIE`/replace_header_files.sed $f > $f~
    mv $f~ $f

done





