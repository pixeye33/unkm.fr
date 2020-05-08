
// first load bounding boxes
var bboxes = null;

var deptData = {};

var wantedCoords = null;

var toBeProcessed = null;

var curShape = null;
var curUnion = null;

var shadow = null;

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

var shadowStyle = {
    "color": "#ffffff",
    "stroke": false,
    "fillOpacity": 0.4
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
    

    // get the shape of the department
    curShape = L.geoJSON(deptData[d], {style: deptStyle});
    curShape.addTo(map);

    // compute a polygonal version of the circle
    var circle = polyCircle([wantedCoords.lng, wantedCoords.lat], 100000.0, 1024);

    // compute union between these two shapes
    var union = martinez.union(deptData[d].geometry.coordinates, [circle]);
    
    // add the union to the rendering
    curUnion = L.geoJSON({
        "type": "Feature",
        "properties": {},
        "geometry": {
            "type": "MultiPolygon",
            "coordinates": union
        }
    }, { style: unionStyle});
    curUnion.addTo(map);
    
    
    var truc = [[[[-90, -360], [90, -360], [90, 360], [-90, 360], [-90, -360]], union[0][0]]].concat(union.slice(1));
    console.log("trcu", JSON.stringify(truc));
    // then draw a shadow outside of the allowed region
    shadow =  L.geoJSON({
        "type": "Feature",
        "properties": {},
        "geometry": {
            "type": "MultiPolygon",
            "coordinates": truc}}, { style: shadowStyle });
    shadow.addTo(map);         
    
    
    // hide the initial big circle
    if (clickBigcircle != undefined) { map.removeLayer(clickBigcircle); };
    
}


function insideShape(coords, shape) {

    if (shape.length == 0)
        return false;
    // a shape is defined by a contour, and possibly holes

    // we first check if the point is in the contour
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

// check if the coordinate is in the shape of the department
function isInDepartment(coords, dept) {
    var polyCoords = deptData[dept]["geometry"]["coordinates"];
    var c = [coords.lng, coords.lat];
    // a shape can be a single shape, or multiple shapes (islands)
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
    
    // if the list of possible shapes is unknown, first load it
    // using the bound boxes and the desired coordinate
    if (toBeProcessed == null) {
        toBeProcessed = getIntersectingBBoxes(wantedCoords);
    }
    
    // then for each possible match
    for(var i = 0; i != toBeProcessed.length; ++i) {
        var d = toBeProcessed[i]["dept"];
        // if the shape is already in memory
        if (d in deptData) {
            // check for intersection
            if (isInDepartment(wantedCoords, d)) {
                // if it is valid, draw this department
                drawDepartment(d);
                break;
            }
        }
        else {
             // otherwise load the data, then restart the display from this department
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
    // clear previous rendering
    if (curShape != null) {
        map.removeLayer(curShape);
    }
    if (curUnion != null) {
        map.removeLayer(curUnion);
    }
    if (shadow != null) {
        map.removeLayer(shadow);
    }
    
    // set the new coordinate
    wantedCoords = e.latlng;
    // if bounding boxes of the departments are ready, 
    // display the department arround the coordinates
    if (bboxes != null) {
        displayDepartmentInternal();
    }
    
}


