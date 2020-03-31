function getNodesInside(elements, center) {
    var result = {};

    for(var i = 0; i != elements.length; ++i) {
        if (elements[i].type == "node") {
            var p = L.latLng(elements[i].lat, elements[i].lon);
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
            if (bestNode < 0 || mAngle < angle) {
                bestNode = n.node;
                angle = mAngle;
            }
        }
    }
    if (bestNode == -1)
        console.log("ERROR: cannot find next point");
    return bestNode;

}


function getGraphLoopFromEdge(line, graph) {
    var result = [];

    var pred = line[0];
    var current = line[1];

    // build the edges
    do {
        var next = findNextOnRegion(pred, current, graph.nodes);
        result.push([current, next]);
        pred = current;
        current = next;
        
    } while ((pred != line[0]) || (current != line[1]));

    return filterBySelectedEdges(graph, result);
}


function getLongitudeCrossing(lat, edge, nodes) {
    var lat1 = nodes[edge[0]].lat;
    var lat2 = nodes[edge[1]].lat;
    var lng1 = nodes[edge[0]].lng;
    var lng2 = nodes[edge[1]].lng;

    if (lat1 > lat2) {
        console.log("ERROR: wrong orientation");
    }

    return lng1 + (lng2 -lng1) / (lat2 - lat1) * (lat - lat1);

    }

    function getContourEdge(graph) {
    var edges = graph.edges;
    var nodes = graph.nodes;
    var lat = window.center.lat;
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

    var contourEdge = getContourEdge(graph);
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
        var edge = graph.edges[e];
        if (edge[0] < edge[1]) {
            var v = edge[0];
            edge[0] = edge[1];
            edge[1] = v;
        }
        var str = edge[0] + "," + edge[1];
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

function setInformationDistance(size) {
    var msg;
    if (size < 1000)
        msg = "Longueur de la boucle&nbsp;: " + Math.round(size) + " m";
    else {
        msg = "Longueur de la boucle&nbsp;: " + (size/1000).toFixed(3).replace(".", ",") + " km";
    }
    document.getElementById("msgLoop").innerHTML = msg;
    setVisibleInformation();
}


function displayCircularPath(e) {
    setVisibleComputation();
    
    removeOldLoop();
   
    var square = clickCircle.getBounds();
    var bboxString = square.getSouth() + "," + square.getWest() + ", " + square.getNorth() + "," + square.getEast();

    var query = "[out:json][timeout:25];(way[\"highway\"][\"highway\"!~\"motorway\"](" + bboxString + ");way[\"pathway\"](" + bboxString + "););out;>;out skel qt;";

    var server = "https://overpass.kumi.systems/api/interpreter";

    var request = new XMLHttpRequest();
    
    
    request.open('GET', server + "?data=" + encodeURIComponent(query), true);
    request.setRequestHeader('Content-type', 'application/json');

    request.addEventListener("load", function() {

        var idProcess = window.idProcess;

        var response = JSON.parse(request.response);

        // get only nodes inside the disc
        var nodes = getNodesInside(response.elements);

        // get all edges inside the disc
        var graph = buildGraphInside(response.elements, nodes);
                                                                
        // select the largest connected component (not perfect to select the main path, but should be ok)
        graph = keepMainCC(graph);

        // build a polygon corresponding to the contour
        graph = getContour(graph);

        if (graph != null) {

            // remove edges twice in the shape
            graph = removeDoubleEdges(graph);

            // select the largest connected component (not perfect to select the main path, but should be ok)
            var mainCC = keepMainCC(graph, true);

            if (mainCC != null) {
                graph = mainCC.graph;
                // draw the contour
                for(var e = 0; e != graph.edges.length; ++e) {
                    window.ClickCircularPaths.push(L.polyline(toPolyLine(graph.edges[e], graph.nodes), { color: '#0060f0'}).addTo(map));
                }

                setInformationDistance(mainCC.size);
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
    if (document.getElementById("loop-panel").style.display == "block") {
        displayCircularPath(e);
    }
    else {
        setVisibleInitialQuestion();
    }
}
