#!/bin/sh

LANG=C

first=1
echo -n "["

for json in [0-9]*.json; do

    if [ $first -eq 1 ]; then
        first=0
    else
        echo -n ","
    fi
    echo ""

    x=$(jq -c .geometry.coordinates $json| tr "\[" " "| tr "\]" "\n"|tr "," " "|grep [0-9]| while read x y; do echo $x; done|datamash min 1 max 1 mean 1)
    y=$(jq -c .geometry.coordinates $json| tr "\[" " "| tr "\]" "\n"|tr "," " "|grep [0-9]| while read x y; do echo $y; done|datamash min 1 max 1 mean 1)
    minx=$(echo $x| cut -d " " -f 1)
    maxx=$(echo $x| cut -d " " -f 2)
    midx=$(echo $x| cut -d " " -f 3)
    miny=$(echo $y| cut -d " " -f 1)
    maxy=$(echo $y| cut -d " " -f 2)
    midy=$(echo $y| cut -d " " -f 3)
    department=$(echo $json| cut -d "." -f1)
    echo -n "\t{ code_insee: \"$department\", minx: \"$minx\", maxx: \"$maxx\", miny: \"$miny\", maxy: \"$maxy\", midx: \"$midx\", midy: \"$midy\"}"
    
done

echo "\n]"
