function getNodesInside(elements, center) {
    var result = {};

    for(var i = 0; i != elements.length; ++i) {
        if (elements[i].type == "node") {
            var p = new L.LatLng(elements[i].lat, elements[i].lon);
            if (p.distanceTo(window.center) < 1000) {
                result[elements[i].id] = p;
                result[elements[i].id].neighbours = [];
                result[elements[i].id].id = elements[i].id;
            }
        }
    }

    return result;
}

function buildGraphInside(elements, nodes) {
    var result = [];
    for(var i = 0; i != elements.length; ++i) {
        if (elements[i].type == "way" || elements[i].type == "pathway") {                
            if (elements[i].nodes.length > 1) {
                for(var n = 1; n != elements[i].nodes.length; ++n) {
                    var idN = elements[i].nodes[n];
                    var idNm1 = elements[i].nodes[n - 1];
                    if (idNm1 in nodes && idN in nodes) {
                        result.push([parseInt(idNm1), parseInt(idN)]);
                        var idE = result.length - 1;
                        nodes[idNm1].neighbours.push({ node: idN, edge: idE});
                        nodes[idN].neighbours.push({ node: idNm1, edge: idE});
                    }
                }
            }
        }
    }

    return {nodes: nodes, edges: result};
}


function filterBySelectedEdges(graph, edges) {
    var newNodes = {};

    // create the set of new nodes
    for(var e = 0; e != edges.length; ++e) {
        newNodes[edges[e][0]] = graph.nodes[edges[e][0]];
        newNodes[edges[e][1]] = graph.nodes[edges[e][1]];
        newNodes[edges[e][0]].neighbours = [];
        newNodes[edges[e][1]].neighbours = [];
    }

    // for each node, create the neighbours list
    for(var e = 0; e != edges.length; ++e) {
        newNodes[edges[e][0]].neighbours.push({ node: edges[e][1], edge: e});
        newNodes[edges[e][1]].neighbours.push({ node: edges[e][0], edge: e});
    }

    return { nodes: newNodes, edges: edges };
}

// compute a pseudo angle defined by a polyline (n1, n2, n3)
function computeAngle(n1, n2, n3) {      
    var v1 = [n2.lat - n1.lat, n2.lng - n1.lng];
    var v2 = [n3.lat - n2.lat, n3.lng - n2.lng];
    var dot = v1[0] * v2[0] + v1[1] * v2[1];
    var det = v1[0] * v2[1] - v1[1] * v2[0];
    return Math.atan2(det, dot);

}

function findNextOnRegion(pred, current, nodes) {
    var predNode = nodes[pred];
    var curNode = nodes[current];
    // if it is a dead-end, we go back to the 
    if (curNode.neighbours.length == 1) {
        return pred;
    }

    var bestNode = -1;
    var angle = 0;
    for(var v = 0; v != curNode.neighbours.length; ++v) {
        var n = curNode.neighbours[v];
        if (pred != n.node) {
            var mAngle = computeAngle(predNode, curNode, nodes[n.node]);
            if (bestNode < 0 || (mAngle < angle)) {
                bestNode = n.node;
                angle = mAngle;
            }
        }
    }
    if (bestNode == -1)
        console.log("ERROR: cannot find next point");
    return bestNode;

}


function intersection(edge1, edge2, graph) {
    
    if ((edge1[0] != edge2[0]) && (edge1[1] != edge2[1]) &&
        (edge1[1] != edge2[0]) && (edge1[0] != edge2[1])) {
        // Inspired from Turfjs packages/turf-line-intersect/index.ts
        const x1 = graph.nodes[edge2[0]].lat;
        const  y1 = graph.nodes[edge2[0]].lng;
        const  x2 = graph.nodes[edge2[1]].lat;
        const  y2 = graph.nodes[edge2[1]].lng;
        const  x3 = graph.nodes[edge1[0]].lat;
        const  y3 = graph.nodes[edge1[0]].lng;
        const  x4 = graph.nodes[edge1[1]].lat;
        const  y4 = graph.nodes[edge1[1]].lng;
        const  denom = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1));
        const  numeA = ((x4 - x3) * (y1 - y3)) - ((y4 - y3) * (x1 - x3));
        const  numeB = ((x2 - x1) * (y1 - y3)) - ((y2 - y1) * (x1 - x3));

        if (denom === 0) {
            if (numeA === 0 && numeB === 0) {
                return null;
            }
            return null;
        }

        const uA = numeA / denom;
        const uB = numeB / denom;

        if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {
            const x = x1 + (uA * (x2 - x1));
            const y = y1 + (uA * (y2 - y1));
            var res = new L.LatLng(x, y);
            res.uA = uA;
            res.uB = uB;
            return res;
        }

    }
    
    return null;
}

function getGraphLoopFromEdge(line, graph) {
    var result = [];

    var pred = line[0];
    var current = line[1];

    // build the edges
    do {
        var next = findNextOnRegion(pred, current, graph.nodes);
        if (next == -1)
            return null;
        result.push([current, next]);

        pred = current;
        current = next;
        
    } while ((pred != line[0]) || (current != line[1]));

    return filterBySelectedEdges(graph, result);
}


function getLongitudeCrossing(lat, edge, nodes) {
    const lat1 = nodes[edge[0]].lat;
    const lat2 = nodes[edge[1]].lat;
    const lng1 = nodes[edge[0]].lng;
    const lng2 = nodes[edge[1]].lng;

    if (lat1 > lat2) {
        console.log("ERROR: wrong orientation");
    }

    return lng1 + (lng2 -lng1) / (lat2 - lat1) * (lat - lat1);

}

function getContourEdge(graph) {
    const edges = graph.edges;
    const nodes = graph.nodes;
    const lat = window.center.lat;
    var crossings = [];

    // get all the edges that are crossing the line defined by the lattitude of the center
    for(var l = 0; l != edges.length; ++l) {
        if (nodes[edges[l][0]].lat <= lat && nodes[edges[l][1]].lat > lat) {
            crossings.push({edge: edges[l], coord: getLongitudeCrossing(lat, edges[l], nodes) });
        }
        else if (nodes[edges[l][1]].lat <= lat && nodes[edges[l][0]].lat > lat) {
            crossings.push({edge: [edges[l][1], edges[l][0]], coord: getLongitudeCrossing(lat, [edges[l][1], edges[l][0]], nodes) });
        }
    }

    if (crossings.length == 0) {
        console.log("ERROR: cannot find a contour edge");
        return null;
    }
    // find the smallest longitude
    var result = 0;
    var longitude = crossings[result].coord;
    for(var i = 0; i < crossings.length; ++i) {
        if (crossings[i].coord < longitude) {
            result = i;
            longitude = crossings[i].coord;
        }
    }

    return crossings[result].edge;

}


function getContour(graph) {

    const contourEdge = getContourEdge(graph);
    if (contourEdge == null)
        return null;

    return getGraphLoopFromEdge(contourEdge, graph);
}

function browseCC(nodeID, graph) {
    graph.nodes[nodeID].seen = true;
    var result = { edges: [], size: 0 };
    for(var nb = 0; nb != graph.nodes[nodeID].neighbours.length; ++nb) {
        var neighbour = graph.nodes[nodeID].neighbours[nb].node;
        var edge = graph.nodes[nodeID].neighbours[nb].edge;
        // to avoid multiple instances of an edge, only add it from the first node
        if (graph.edges[edge][0] == nodeID)
            result.edges.push(graph.edges[edge]);
        // then use recursive call if required
        if (!graph.nodes[neighbour].seen) {  
            result.size += graph.nodes[nodeID].distanceTo(graph.nodes[neighbour]);
            var run = browseCC(neighbour, graph);
            result.edges = result.edges.concat(run.edges);
            result.size += run.size;
        }
    }

    return result;
}

function keepMainCC(graph, withSize = false) {
    var cc = [];

    for(var n in graph.nodes) {
        graph.nodes[n].seen = false;
    }


    // compute connecte component size
    for(var n in graph.nodes) {
        n = parseInt(n);
        if (!graph.nodes[n].seen) {
            cc.push(browseCC(n, graph));
        }
    }

    if (cc.length == 0)
        return null;

    // select the largest one
    var idCC = 0;
    var size = cc[idCC].size;
    for(var i = 1; i != cc.length; ++i) {
        if (cc[i].size > size) {
            idCC = i;
            size = cc[i].size;
        }
    }
    var result = filterBySelectedEdges(graph, cc[idCC].edges);

    if (withSize) {
        return { graph: result, size: size };
    }
    else
        return result;
}


function removeDoubleEdges(graph) {
    var edges = {};

    for(var e = 0; e != graph.edges.length; ++e) {
        const edge = graph.edges[e];
        if (edge[0] < edge[1]) {
            var v = edge[0];
            edge[0] = edge[1];
            edge[1] = v;
        }
        const str = edge[0] + "," + edge[1];
        if (str in edges) {
            delete edges[str];
        }
        else {
            edges[str] = edge;
        }
    }
    var flatEdges = [];
    for(var e in edges) {
        flatEdges.push(edges[e]);
    }

    return filterBySelectedEdges(graph, flatEdges);
}

function toPolyLine(edge, nodes) {
    return [[nodes[edge[0]].lat, nodes[edge[0]].lng], [nodes[edge[1]].lat, nodes[edge[1]].lng]];
}


function removeOldLoop() {
     if (window.ClickCircularPaths.length != 0) { 
        for(var p = 0; p != window.ClickCircularPaths.length; ++p) {
            map.removeLayer(window.ClickCircularPaths[p]); 
        }
    }
}

function setInformationError() {
    document.getElementById("msgLoop").innerHTML = "Erreur de calcul... Désolé :(";
    setVisibleInformation();
}

function distToStr(dist) {
    if (dist < 1000)
        return Math.round(dist) + " m";
    else {
        return  (dist/1000).toFixed(3).replace(".", ",") + " km";
    }
}
function setInformationDistance(size, sizePath) {
    var msg = "Longueur de la boucle&nbsp;: <strong>" + distToStr(size) + "</strong>";
    if (sizePath > 0) {
        msg += ".<br /> Longueur du chemin pour s'y rendre&nbsp;: <strong>" + distToStr(sizePath) + ".</strong>";
    }
    
    document.getElementById("msgLoop").innerHTML = msg;
    setVisibleInformation();
}

function newIDPoint(graph, startingFrom = 0) {
    var i = startingFrom;
    while(i in graph.nodes) {
        ++i;
    }
    return i;
}

function mergeSameLocation(graph) {
    var toBeDeleted = [];
    
    var list = {};
    for(var n1 in graph.nodes) {
        const str = graph.nodes[n1].lat + "," + graph.nodes[n1].lng;
        if (str in list) {
            toBeDeleted.push(n1);
        }
    }
    
    for(var i = 0; i != toBeDeleted.length; ++i)
        delete graph.nodes[i];
        
    // clear nodes
    for(var e = 0; e != graph.edges.length; ++e) {
        graph.nodes[graph.edges[e][0]].neighbours = [];
        graph.nodes[graph.edges[e][1]].neighbours = [];
    }

    // for each node, create the neighbours list
    for(var e = 0; e != graph.edges.length; ++e) {
        graph.nodes[graph.edges[e][0]].neighbours.push({ node: graph.edges[e][1], edge: e});
        graph.nodes[graph.edges[e][1]].neighbours.push({ node: graph.edges[e][0], edge: e});
    }

    return graph;    
}

function boxFromEdge(edge, graph) {
    const p1 = graph.nodes[edge[0]];
    const p2 = graph.nodes[edge[1]];
    const minX = p1.lat < p2.lat ? p1.lat : p2.lat;
    const minY = p1.lng < p2.lng ? p1.lng : p2.lng;
    const maxX = p1.lat > p2.lat ? p1.lat : p2.lat;
    const maxY = p1.lng > p2.lng ? p1.lng : p2.lng;
    return {minX: minX, minY: minY, maxX: maxX, maxY: maxY};
}

function addMiddleCrossing(graph) {
     
    const index = new Flatbush(graph.edges.length);
    for (const e of graph.edges) {
        const b = boxFromEdge(e, graph);
        index.add(b.minX, b.minY, b.maxX, b.maxY);
    }
    index.finish();
    
    // compute intersections
    
    var intersections = [];
    for(var e1 = 0; e1 != graph.edges.length; ++e1) {
        const b = boxFromEdge(graph.edges[e1], graph);
        const found = index.search(b.minX, b.minY, b.maxX, b.maxY);
        for(var e2 of found) {
            if (e1 > e2) {
                const inter = intersection(graph.edges[e1], graph.edges[e2], graph);
                if (inter != null) {
                    intersections.push({e1: e1, e2: e2, coord: inter});
                }
            }
        }
    }
    
    // add nodes corresponding to these intersections and group them in each old edge
    var edgeSplit = {};
    for(const inter of intersections) {
        if (inter !== undefined) {
            var newID = newIDPoint(graph, inter.e1);
            graph.nodes[newID] = new L.LatLng(inter.coord.lat, inter.coord.lng);
            graph.nodes[newID].id = newID;
            if (!(inter.e1 in edgeSplit)) edgeSplit[inter.e1] = [];
            if (!(inter.e2 in edgeSplit)) edgeSplit[inter.e2] = [];
            edgeSplit[inter.e1].push({id: newID, coord: inter.coord.uA, origin: [e1, e2]});
            edgeSplit[inter.e2].push({id: newID, coord: inter.coord.uB, origin: [e1, e2]});
        }
    }
    
    // add missing edges along each initial edge
    for(var es in edgeSplit) {
        if (edgeSplit[es].length != 0) {
            const prev = graph.edges[es][1];
            graph.edges[es][1] = edgeSplit[es][0].id;            
            if (edgeSplit[es].length > 1) {
                edgeSplit[es].sort(function(a, b) { return a.coord - b.coord; });
                for(var e = 1; e != edgeSplit[es].length; ++e) {
                    var p1 = edgeSplit[es][e - 1].id;
                    var p2 = edgeSplit[es][e].id;
                    graph.edges.push([p1, p2]);
                }
            }
            graph.edges.push([edgeSplit[es][edgeSplit[es].length - 1].id, prev]);
        }
    }

    
    // clear nodes
    for(var e = 0; e != graph.edges.length; ++e) {
        graph.nodes[graph.edges[e][0]].neighbours = [];
        graph.nodes[graph.edges[e][1]].neighbours = [];
    }

    // for each node, create the neighbours list
    for(var e = 0; e != graph.edges.length; ++e) {
        graph.nodes[graph.edges[e][0]].neighbours.push({ node: graph.edges[e][1], edge: e});
        graph.nodes[graph.edges[e][1]].neighbours.push({ node: graph.edges[e][0], edge: e});
    }

    return graph;
   
}

function addAutointersections(graph) {
    graph = mergeSameLocation(graph);
    
    graph = addMiddleCrossing(graph);
    
    return graph;
}

// return distance between a point C and a segment [A, B]
function distancePointEdge(C, A, B) {
    // cf http://www.faqs.org/faqs/graphics/algorithms-faq/
    // Subject 1.02: How do I find the distance from a point to a line?
    
    const dist = B.distanceTo(A);
    if (dist == 0)
        return {distance: C.distanceTo(A), point: A};
    const dist2 = dist * dist;
    
    const r = ((C.x - A.x) * (B.x - A.x) +
               (C.y - A.y) * (B.y - A.y)) / dist2;

    if (r < 0) {
        return {distance: C.distanceTo(A), point: A};
    }
    else if (r > 1) {
        return {distance: C.distanceTo(B), point: B};
    }
    else {
        const Px = A.x + r * (B.x - A.x);
        const Py = A.y + r * (B.y - A.y);
        const point = new L.Point(Px, Py);
        
        return {distance: C.distanceTo(point), point: point};
    }
}


// get the point of graph closest to given point. 
// the returned point can be an existing node, or a coordinate along an edge
function getClosestPoint(graph, point) {
    var dist = -1;
    var cPoint;
    var eID;

    var pPoint = L.Projection.SphericalMercator.project(point);
    for(var e = 0; e != graph.edges.length; ++e) {
        var A = L.Projection.SphericalMercator.project(graph.nodes[graph.edges[e][0]]);
        var B = L.Projection.SphericalMercator.project(graph.nodes[graph.edges[e][1]]);
        var cp = distancePointEdge(pPoint, A, B);
        if (dist < 0 || cp.distance < dist) {
            dist = cp.distance;
            cPoint = cp.point;
            eID = e;
        }
    }
    cPoint = L.Projection.SphericalMercator.unproject(cPoint);
    
    return {distance: dist, point: cPoint, edgeID: eID};
    
}

// get the node of graph closest to point, adding it on an edge if required
function getClosestNode(graph, point) {
    var cPoint = getClosestPoint(graph, point);
    
    if ('id' in cPoint.point) {
        return cPoint.point;
    }
    else {
        
        var newID = newIDPoint(graph);
        graph.nodes[newID] = cPoint.point;
        graph.nodes[newID].id = newID;
        var edge = graph.edges[cPoint.edgeID];
        graph.edges.push([edge[0], newID]);
        graph.edges[cPoint.edgeID][0] = newID;
        
        // update graph
        graph = filterBySelectedEdges(graph, graph.edges);
        
        return graph.nodes[newID];
    }
}

function addPathToPoint(graph, point) {
    // get the node of graph closest to point, adding it on an edge if required
    var cNode = getClosestNode(graph, point);
 
    // if this point is exactly the request point, we return it
    if (cNode.distanceTo(point) == 0) {
        return { graph: graph, point: cNode.id};
    }
    else {
        // otherwise we add a new one
        var newID = newIDPoint(graph);
        graph.nodes[newID] = new L.LatLng(point.lat, point.lng);
        graph.nodes[newID].id = newID;
        graph.edges.push([newID, cNode.id]);
        graph = filterBySelectedEdges(graph, graph.edges);
        return { graph: graph, point: newID};
    }
}

// Dijkstra implementation

function findShortestPath(graph, source, targets) {
    // imported from: https://github.com/mburst/dijkstras-algorithm/blob/master/dijkstras.js
    function PriorityQueue () {
        this._nodes = [];

        this.enqueue = function (priority, key) {
            this._nodes.push({key: key, priority: priority });
            this.sort();
        };
        this.dequeue = function () {
            return this._nodes.shift().key;
        };
        this.sort = function () {
            this._nodes.sort(function (a, b) {
            return a.priority - b.priority;
            });
        };
        this.isEmpty = function () {
            return !this._nodes.length;
        };
    }
    const INFINITY = 1/0;


    
    var path = [];
    
    if (targets.indexOf(source) != -1)
        return path;
    
    var nodes = new PriorityQueue(),
        distances = {},
        previous = {},
        smallest, vertex, neighbor, alt;

    for(var n in graph.nodes) {
      if(n == source) {
        distances[n] = 0;
        nodes.enqueue(0, n);
      }
      else {
        distances[n] = INFINITY;
        nodes.enqueue(INFINITY, n);
      }

      previous[n] = null;
    }

    while(!nodes.isEmpty()) {
      smallest = nodes.dequeue();
      

      if(targets.indexOf(smallest) != -1) {
        path = [];

        while(previous[smallest]) {
          path.push(smallest);
          smallest = previous[smallest];
        }
        path.push(source);

        break;
      }

      if(!smallest || distances[smallest] === INFINITY){
        continue;
      }

      for(var n = 0; n != graph.nodes[smallest].neighbours.length; ++n) {
        var neighbor = graph.nodes[smallest].neighbours[n];
        var nodeID = "" + neighbor.node;
        alt = distances[smallest] + graph.nodes[smallest].distanceTo(graph.nodes[nodeID]);

        if(alt < distances[nodeID]) {
          distances[nodeID] = alt;
          previous[nodeID] = smallest;
          nodes.enqueue(alt, nodeID);
        }
      }
    }
    return path;
    
}

function computeShortestPath(graph, targetNodes) {
    // find the closest point between the selected point and the graph
    // and add the corresponding path
    var size = graph.edges.length;
    var pointAndGraph = addPathToPoint(graph, window.center);
    var graph = pointAndGraph.graph;
    var source = pointAndGraph.point;
    
    // find the shortest path between source and targetNodes in graph
    var targets = [];
    for(var n in targetNodes) {
        targets.push(n);
    }
    
    
    // if a supplementary point has been added on the graph
    if (graph.edges.length != size) {
        // and if this new point is not the center
        if (graph.nodes[source].neighbours.length == 1) {
            // add it to the target list
            targets.push(graph.nodes[source].neighbours[0].node);
        }
    }
    
    
    // compute shortest path between source and targets
    var path = findShortestPath(graph, source, targets);
    
    // build list of edges
    var edges = [];
    if (path.length > 1) {
        for(var p = 1; p != path.length; ++p) {
            edges.push([path[p - 1], path[p]]);
        }
    }
    
    // return the corresponding graph
    return filterBySelectedEdges(graph, edges);
}

function copyGraph(graph) {
    var result = {
       edges: [],
       nodes: {}
    };
    
    result.edges = JSON.parse(JSON.stringify(graph.edges));
    
    for(var n in graph.nodes) {
        result.nodes[n] = new L.LatLng(graph.nodes[n].lat, graph.nodes[n].lng);
        result.nodes[n].id = graph.nodes[n].id;
        result.nodes[n].neighbours = [];
        for(var nb = 0; nb != graph.nodes[n].neighbours.length; ++nb) {
            result.nodes[n].neighbours.push(JSON.parse(JSON.stringify(graph.nodes[n].neighbours[nb])));
        }
    }
    return result;
    
}

function computeSize(graph) {
    var result = 0;
    
    for(var edge of graph.edges) {
        result += graph.nodes[edge[0]].distanceTo(graph.nodes[edge[1]]);
    }
    
    return result;
    
}

function displayCircularPath(e) {
    setVisibleComputation();
    
    removeOldLoop();
   
    const square = clickCircle.getBounds();
    const bboxString = square.getSouth() + "," + square.getWest() + ", " + square.getNorth() + "," + square.getEast();

    const query = "[out:json][timeout:25];(way[\"highway\"][\"highway\"!~\"motorway\"](" + bboxString + ");way[\"pathway\"](" + bboxString + "););out;>;out skel qt;";

    const server = "https://overpass.kumi.systems/api/interpreter";

    const request = new XMLHttpRequest();
    
    
    request.open('GET', server + "?data=" + encodeURIComponent(query), true);
    request.setRequestHeader('Content-type', 'application/json');

    request.addEventListener("load", function() {

        const idProcess = window.idProcess;

        const response = JSON.parse(request.response);


        // get only nodes inside the disc
        var nodes = getNodesInside(response.elements);

        // get all edges inside the disc
        var graph = buildGraphInside(response.elements, nodes);
        
        // create a copy of the graph to find shortest path
        var fullGraph = copyGraph(graph);
                                                                
        // select the largest connected component (not perfect to select the main path, but should be ok)
        graph = keepMainCC(graph);
        
        // add self intersections as new points: it is a trick to avoid 
        // inversions, that is not satisfying (for example when a self intersection
        // corresponds to a bridge), but cannot handle contour without that trick (what is inside 
        // and outside?)
        graph = addAutointersections(graph);

        // build a polygon corresponding to the contour
        graph = getContour(graph);

        if (graph != null) {

            // remove edges twice in the shape
            graph = removeDoubleEdges(graph);

            // select the largest connected component (not perfect to select the main path, but should be ok)
            const mainCC = keepMainCC(graph, true);

            if (mainCC != null) {
                graph = mainCC.graph;
                
                // build a copy of the shortest path
                var shortestPath = computeShortestPath(fullGraph, graph.nodes);
                


                // draw the shortest path
                for(var e = 0; e != shortestPath.edges.length; ++e) {
                    window.ClickCircularPaths.push(L.polyline(toPolyLine(shortestPath.edges[e], shortestPath.nodes), { color: '#2ab50a'}).addTo(map));
                }

                // draw the contour
                for(var e = 0; e != graph.edges.length; ++e) {
                    window.ClickCircularPaths.push(L.polyline(toPolyLine(graph.edges[e], graph.nodes), { color: '#0060f0'}).addTo(map));
                }
                
                setInformationDistance(mainCC.size, computeSize(shortestPath));
            }
            else {
                setInformationError();
            }
        }
        else {
            setInformationError();
        }
        
        

    });
    request.send();

}

function locationAndLoopUpdate(e) {
    if (document.getElementById("loop").className.split(" ").indexOf("openPanel") != -1) {
        if (typeof window.location.hash != "undefined" && window.location.hash.length > 5) {
            displayCircularPath(e);
        }
        else {
            setVisibleNoDisc();
        }
            
    }
    else {
        setVisibleInitialQuestion();
    }
}
