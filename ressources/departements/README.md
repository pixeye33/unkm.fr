# Description des données

All data contained in this directory has been transformed from the administratives boundaries of the french departements in OpenStreetMap, thus are distributed under [ODbL license](https://www.openstreetmap.org/copyright).
They were [extracted](http://osm13.openstreetmap.fr/~cquest/openfla/export/) on 11 Fabruary 2019.

These geometric shapes are too complex, thus they were simplified as suggested by  [Grégoire David](https://github.com/gregoiredavid/france-geojson) with the following command line:

    mapshaper -i input.shp snap -proj wgs84 -simplify 25% weighted keep-shapes -o format=geojson precision=0.00001 output.json

This first file has been split into one file per department, in order to reduce the data to be downloaded by each unkm.fr user. We used the following command line:

    jq -c ".features[] | .properties.code_insee" output.json | while read line; do fname=${line//\"}.json; jq -c ".features[] | select(.properties.code_insee==$line)" output.json > $fname; done


Finally, a small script (bboxes.sh) can be run to produce a file with the bounding boxes of each department.
