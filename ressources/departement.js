
// first load bounding boxes
var bboxes = null;

var deptData = {};

var wantedCoords = null;

var toBeProcessed = null;

var curShape = null;
var curUnion = null;


var deptStyle = {
    "color": "#EC000C",
    "weight": 3,
    "opacity": 0.5,
    "fill": false,
    "dashArray": '6, 6'
};

var unionStyle = {
    "color": "#EC000C",
    "weight": 3,
    "opacity": 1,
    "fill": false
};

function loadJSON(callback) {   

    var xobj = new XMLHttpRequest();
        xobj.overrideMimeType("application/json");
    xobj.open('GET', 'ressources/departements/bboxes.json', true);
    xobj.onreadystatechange = function () {
          if (xobj.readyState == 4 && xobj.status == "200") {
            // Required use of an anonymous callback as .open will NOT return a value but simply returns undefined in asynchronous mode
            callback(xobj.responseText);
          }
    };
    xobj.send(null);  
}

loadJSON(function(response) {
   bboxes = JSON.parse(response);
   if (wantedCoords != null) {
        displayDepartmentInternal();
   }
});

function distance2(p1, p2) {
    return (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]);
}

function getIntersectingBBoxes(latlng) {
    var result = [];
    for(var ll of bboxes["bboxes"]) {
        if (ll["miny"] <= latlng.lat && ll["minx"] <= latlng.lng &&
            ll["maxy"] >= latlng.lat && ll["maxx"] >= latlng.lng) {
            result.push({ "dept": ll["code_insee"], "distance":
                distance2([latlng.lng, latlng.lat], [ll["midx"], ll["midy"]])});
        }
    }
    
    
    result = result.sort(function(a, b) { return a["distance"] - b["distance"]; });
    
    return result;
}

function inside(point, vs) {
    // function from  https://github.com/substack/point-in-polygon (MIT license)
    // ray-casting algorithm based on
    // http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

    var x = point[0], y = point[1];
    
    var inside = false;
    for (var i = 0, j = vs.length - 1; i < vs.length; j = i++) {
        var xi = vs[i][0], yi = vs[i][1];
        var xj = vs[j][0], yj = vs[j][1];

        var intersect = ((yi > y) != (yj > y))
            && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }
    
    return inside;
};


function loadDepartment(d, callback) {
    var xobj = new XMLHttpRequest();
        xobj.overrideMimeType("application/json");
    xobj.open('GET', 'ressources/departements/' + d + '.json', true);
    xobj.onreadystatechange = function () {
          if (xobj.readyState == 4 && xobj.status == "200") {
            // Required use of an anonymous callback as .open will NOT return a value but simply returns undefined in asynchronous mode
            callback(xobj.responseText);
          }
    };
    xobj.send(null);  
}


// circle discretization from
// https://github.com/gabzim/circle-to-polygon/blob/master/index.js
function toRadians(angleInDegrees) {
  return (angleInDegrees * Math.PI) / 180;
}

function toDegrees(angleInRadians) {
  return (angleInRadians * 180) / Math.PI;
}

function offset(c1, distance, bearing) {
  var lat1 = toRadians(c1[1]);
  var lon1 = toRadians(c1[0]);
  
  var dByR = distance / 6378137; // distance divided by 6378137 (radius of the earth) wgs84
  
  console.log("dByR 3", dByR);
  
  var lat = Math.asin(
    Math.sin(lat1) * Math.cos(dByR) +
      Math.cos(lat1) * Math.sin(dByR) * Math.cos(bearing)
  );
  var lon =
    lon1 +
    Math.atan2(
      Math.sin(bearing) * Math.sin(dByR) * Math.cos(lat1),
      Math.cos(dByR) - Math.sin(lat1) * Math.sin(lat)
    );
  return [toDegrees(lon), toDegrees(lat)];
}

function polyCircle(center, radius, nbSegments) {
    var n = nbSegments ? nbSegments : 32;
    var coordinates = [];
    for (var i = 0; i < n; ++i) {
        coordinates.push(offset(center, radius, (2 * Math.PI * -i) / n));
    }
    coordinates.push(coordinates[0]);   
    
    return coordinates;
}

function drawDepartment(d) {
    
    // hide the initial big circle
    if (clickBigcircle != undefined) { map.removeLayer(clickBigcircle); };
    
    curShape = L.geoJSON(deptData[d], {style: deptStyle});
    curShape.addTo(map);


    var circle = polyCircle([wantedCoords.lng, wantedCoords.lat], 100000.0, 1024);

    
    var union = martinez.union(deptData[d].geometry.coordinates, [circle]);
    console.log("union", JSON.stringify(union));
    
    curUnion = L.geoJSON({
        "type": "Feature",
        "properties": {},
        "geometry": {
            "type": "MultiPolygon",
            "coordinates": union
        }
    }, { style: unionStyle});
    curUnion.addTo(map);
    
    
}

function insideShape(coords, shape) {
    if (shape.length == 0)
        return false;
    if (inside(coords, shape[0])) {
        if (shape.length == 1) {
            // no hole: we are "inside"
            return true;
        }
        else {
            // check the holes
            for(var i = 1; i != shape.length; ++i) {
                if (inside(coords, shape[i]))
                    // in a hole means outside of the shape
                    return false;
            }
            return true;
        }
    }
    else {
        return false;
    }
    
}

function isInDepartment(coords, dept) {
    var polyCoords = deptData[dept]["geometry"]["coordinates"];
    var c = [coords.lng, coords.lat];
    if (deptData[dept]["geometry"]["type"] == "Polygon") {
        return insideShape(c, polyCoords);
    }
    else {
        if (deptData[dept]["geometry"]["type"] == "MultiPolygon") {
            for(var p of polyCoords) {
                if (insideShape(c, p)) {
                    return true;
                }
            }
            return false;
        }
        else {
            return false;
        }
    }
}

function displayDepartmentInternal() {
    
    if (toBeProcessed == null) {
        toBeProcessed = getIntersectingBBoxes(wantedCoords);
    }
    
    for(var i = 0; i != toBeProcessed.length; ++i) {
        var d = toBeProcessed[i]["dept"];
        if (d in deptData) {
            if (isInDepartment(wantedCoords, d)) {
                drawDepartment(d);
                break;
            }
        }
        else {
             toBeProcessed = toBeProcessed.slice(i);
             loadDepartment(d, function(response) {
                 deptData[d] = JSON.parse(response);
                 displayDepartmentInternal();
             });
             return;
        }
    }
    
    wantedCoords = null;
    toBeProcessed = null;
}

function displayDepartment(e) {
    if (curShape != null) {
        map.removeLayer(curShape);
    }
    if (curUnion != null) {
        map.removeLayer(curUnion);
    }
    wantedCoords = e.latlng;
    if (bboxes != null) {
        displayDepartmentInternal();
    }
    
}


