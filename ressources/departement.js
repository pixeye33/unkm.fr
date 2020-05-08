
// first load bounding boxes
var bboxes = null;

var deptData = {};

var wantedCoords = null;

var toBeProcessed = null;

var curShape = null;

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

function drawDepartment(d) {
    curShape = L.geoJSON(deptData[d]);
    curShape.addTo(map);
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
    wantedCoords = e.latlng;
    if (bboxes != null) {
        displayDepartmentInternal();
    }
    
}


