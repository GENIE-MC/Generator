shopt -s nullglob
for f in *.*
do
    echo "file - $f"
    sed -f /hepstore/mroda/Software/Genie/generator/Sources/v3config2/replace_header_files.sed $f > $f~
    mv $f~ $f

done





