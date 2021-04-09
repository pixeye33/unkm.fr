/** This first part comes from flatbush.min.js 
    https://github.com/mourner/flatbush/ 
    a library for fast search in 2D plan
 */

!function(t,i){"object"==typeof exports&&"undefined"!=typeof module?module.exports=i():"function"==typeof define&&define.amd?define(i):(t=t||self).Flatbush=i()}(this,(function(){"use strict";var t=function(){this.ids=[],this.values=[],this.length=0};t.prototype.clear=function(){this.length=0},t.prototype.push=function(t,i){var s=this.length++;for(this.ids[s]=t,this.values[s]=i;s>0;){var e=s-1>>1,h=this.values[e];if(i>=h)break;this.ids[s]=this.ids[e],this.values[s]=h,s=e}this.ids[s]=t,this.values[s]=i},t.prototype.pop=function(){if(0!==this.length){var t=this.ids[0];if(this.length--,this.length>0){for(var i=this.ids[0]=this.ids[this.length],s=this.values[0]=this.values[this.length],e=this.length>>1,h=0;h<e;){var r=1+(h<<1),n=r+1,o=this.ids[r],a=this.values[r],u=this.values[n];if(n<this.length&&u<a&&(r=n,o=this.ids[n],a=u),a>=s)break;this.ids[h]=o,this.values[h]=a,h=r}this.ids[h]=i,this.values[h]=s}return t}},t.prototype.peek=function(){return this.ids[0]},t.prototype.peekValue=function(){return this.values[0]};var i=[Int8Array,Uint8Array,Uint8ClampedArray,Int16Array,Uint16Array,Int32Array,Uint32Array,Float32Array,Float64Array],s=function(s,e,h,r){if(void 0===e&&(e=16),void 0===h&&(h=Float64Array),void 0===s)throw new Error("Missing required argument: numItems.");if(isNaN(s)||s<=0)throw new Error("Unpexpected numItems value: "+s+".");this.numItems=+s,this.nodeSize=Math.min(Math.max(+e,2),65535);var n=s,o=n;this._levelBounds=[4*n];do{o+=n=Math.ceil(n/this.nodeSize),this._levelBounds.push(4*o)}while(1!==n);this.ArrayType=h||Float64Array,this.IndexArrayType=o<16384?Uint16Array:Uint32Array;var a=i.indexOf(this.ArrayType),u=4*o*this.ArrayType.BYTES_PER_ELEMENT;if(a<0)throw new Error("Unexpected typed array class: "+h+".");r&&r instanceof ArrayBuffer?(this.data=r,this._boxes=new this.ArrayType(this.data,8,4*o),this._indices=new this.IndexArrayType(this.data,8+u,o),this._pos=4*o,this.minX=this._boxes[this._pos-4],this.minY=this._boxes[this._pos-3],this.maxX=this._boxes[this._pos-2],this.maxY=this._boxes[this._pos-1]):(this.data=new ArrayBuffer(8+u+o*this.IndexArrayType.BYTES_PER_ELEMENT),this._boxes=new this.ArrayType(this.data,8,4*o),this._indices=new this.IndexArrayType(this.data,8+u,o),this._pos=0,this.minX=1/0,this.minY=1/0,this.maxX=-1/0,this.maxY=-1/0,new Uint8Array(this.data,0,2).set([251,48+a]),new Uint16Array(this.data,2,1)[0]=e,new Uint32Array(this.data,4,1)[0]=s),this._queue=new t};function e(t,i,s){return t<i?i-t:t<=s?0:t-s}function h(t,i){for(var s=0,e=i.length-1;s<e;){var h=s+e>>1;i[h]>t?e=h:s=h+1}return i[s]}function r(t,i,s,e,h){var r=t[e];t[e]=t[h],t[h]=r;var n=4*e,o=4*h,a=i[n],u=i[n+1],d=i[n+2],p=i[n+3];i[n]=i[o],i[n+1]=i[o+1],i[n+2]=i[o+2],i[n+3]=i[o+3],i[o]=a,i[o+1]=u,i[o+2]=d,i[o+3]=p;var f=s[e];s[e]=s[h],s[h]=f}function n(t,i){var s=t^i,e=65535^s,h=65535^(t|i),r=t&(65535^i),n=s|e>>1,o=s>>1^s,a=h>>1^e&r>>1^h,u=s&h>>1^r>>1^r;o=(s=n)&(e=o)>>2^e&(s^e)>>2,a^=s&(h=a)>>2^e&(r=u)>>2,u^=e&h>>2^(s^e)&r>>2,o=(s=n=s&s>>2^e&e>>2)&(e=o)>>4^e&(s^e)>>4,a^=s&(h=a)>>4^e&(r=u)>>4,u^=e&h>>4^(s^e)&r>>4,a^=(s=n=s&s>>4^e&e>>4)&(h=a)>>8^(e=o)&(r=u)>>8;var d=t^i,p=(e=(u^=e&h>>8^(s^e)&r>>8)^u>>1)|65535^(d|(s=a^a>>1));return((p=1431655765&((p=858993459&((p=252645135&((p=16711935&(p|p<<8))|p<<4))|p<<2))|p<<1))<<1|(d=1431655765&((d=858993459&((d=252645135&((d=16711935&(d|d<<8))|d<<4))|d<<2))|d<<1)))>>>0}return s.from=function(t){if(!(t instanceof ArrayBuffer))throw new Error("Data must be an instance of ArrayBuffer.");var e=new Uint8Array(t,0,2),h=e[0],r=e[1];if(251!==h)throw new Error("Data does not appear to be in a Flatbush format.");if(r>>4!=3)throw new Error("Got v"+(r>>4)+" data when expected v3.");var n=new Uint16Array(t,2,1)[0],o=new Uint32Array(t,4,1)[0];return new s(o,n,i[15&r],t)},s.prototype.add=function(t,i,s,e){var h=this._pos>>2;return this._indices[h]=h,this._boxes[this._pos++]=t,this._boxes[this._pos++]=i,this._boxes[this._pos++]=s,this._boxes[this._pos++]=e,t<this.minX&&(this.minX=t),i<this.minY&&(this.minY=i),s>this.maxX&&(this.maxX=s),e>this.maxY&&(this.maxY=e),h},s.prototype.finish=function(){if(this._pos>>2!==this.numItems)throw new Error("Added "+(this._pos>>2)+" items when expected "+this.numItems+".");for(var t=this.maxX-this.minX,i=this.maxY-this.minY,s=new Uint32Array(this.numItems),e=0;e<this.numItems;e++){var h=4*e,o=this._boxes[h++],a=this._boxes[h++],u=this._boxes[h++],d=this._boxes[h++],p=Math.floor(65535*((o+u)/2-this.minX)/t),f=Math.floor(65535*((a+d)/2-this.minY)/i);s[e]=n(p,f)}!function t(i,s,e,h,n){if(!(h>=n)){for(var o=i[h+n>>1],a=h-1,u=n+1;;){do{a++}while(i[a]<o);do{u--}while(i[u]>o);if(a>=u)break;r(i,s,e,a,u)}t(i,s,e,h,u),t(i,s,e,u+1,n)}}(s,this._boxes,this._indices,0,this.numItems-1);for(var _=0,l=0;_<this._levelBounds.length-1;_++)for(var v=this._levelBounds[_];l<v;){for(var x=1/0,y=1/0,m=-1/0,c=-1/0,b=l,w=0;w<this.nodeSize&&l<v;w++){var A=this._boxes[l++],g=this._boxes[l++],E=this._boxes[l++],I=this._boxes[l++];A<x&&(x=A),g<y&&(y=g),E>m&&(m=E),I>c&&(c=I)}this._indices[this._pos>>2]=b,this._boxes[this._pos++]=x,this._boxes[this._pos++]=y,this._boxes[this._pos++]=m,this._boxes[this._pos++]=c}},s.prototype.search=function(t,i,s,e,h){if(this._pos!==this._boxes.length)throw new Error("Data not yet indexed - call index.finish().");for(var r=this._boxes.length-4,n=this._levelBounds.length-1,o=[],a=[];void 0!==r;){for(var u=Math.min(r+4*this.nodeSize,this._levelBounds[n]),d=r;d<u;d+=4){var p=0|this._indices[d>>2];s<this._boxes[d]||e<this._boxes[d+1]||t>this._boxes[d+2]||i>this._boxes[d+3]||(r<4*this.numItems?(void 0===h||h(p))&&a.push(p):(o.push(p),o.push(n-1)))}n=o.pop(),r=o.pop()}return a},s.prototype.neighbors=function(t,i,s,r,n){if(void 0===s&&(s=1/0),void 0===r&&(r=1/0),this._pos!==this._boxes.length)throw new Error("Data not yet indexed - call index.finish().");for(var o=this._boxes.length-4,a=this._queue,u=[],d=r*r;void 0!==o;){for(var p=Math.min(o+4*this.nodeSize,h(o,this._levelBounds)),f=o;f<p;f+=4){var _=0|this._indices[f>>2],l=e(t,this._boxes[f],this._boxes[f+2]),v=e(i,this._boxes[f+1],this._boxes[f+3]),x=l*l+v*v;o<4*this.numItems?(void 0===n||n(_))&&a.push(-_-1,x):a.push(_,x)}for(;a.length&&a.peek()<0;){if(a.peekValue()>d)return a.clear(),u;if(u.push(-a.pop()-1),u.length===s)return a.clear(),u}o=a.pop()}return a.clear(),u},s}));


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


const gpxHeader = `<?xml version="1.0" encoding="UTF-8"?>
<gpx
  version="1.0"
  creator="GPX from IGN benchmarks list"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
  xmlns="http://www.topografix.com/GPX/1/0"
  xmlns:locus="http://www.locusmap.eu"
  xsi:schemaLocation="http://www.topografix.com/GPX/1/0 http://www.topografix.com/GPX/1/0/gpx.xsd">
  <name>Boucle de marche en confinement</name>
  <desc>Plus grande boucle dans un rayon de 1km</desc>
`;

const gpxFooter = `
</gpx>
`;

var gpxBody = "";


/*!
 * JavaScript function to calculate the destination point given start point latitude / longitude (numeric degrees), bearing (numeric degrees) and distance (in m).
 *
 * Taken from http://movable-type.co.uk/scripts/latlong-vincenty-direct.html and optimized / cleaned up by Mathias Bynens <http://mathiasbynens.be/>
 * Source : https://gist.github.com/mathiasbynens/354587/4137c2350d7cee9b757444ef3a6f13a1f69c2abc
 * Based on the Vincenty direct formula by T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975 <http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf>
 */

function toRad(degrees)
{
    return degrees * (Math.PI/180);
}
function toDeg(radians)
{
    return radians * (180 / Math.PI);
}

function destVincenty(lat1, lon1, brng, dist) {
    var a = 6378137,
        b = 6356752.3142,
        f = 1 / 298.257223563, // WGS-84 ellipsiod
        s = dist,
        alpha1 = toRad(brng),
        sinAlpha1 = Math.sin(alpha1),
        cosAlpha1 = Math.cos(alpha1),
        tanU1 = (1 - f) * Math.tan(toRad(lat1)),
        cosU1 = 1 / Math.sqrt((1 + tanU1 * tanU1)), sinU1 = tanU1 * cosU1,
        sigma1 = Math.atan2(tanU1, cosAlpha1),
        sinAlpha = cosU1 * sinAlpha1,
        cosSqAlpha = 1 - sinAlpha * sinAlpha,
        uSq = cosSqAlpha * (a * a - b * b) / (b * b),
        A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq))),
        B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq))),
        sigma = s / (b * A),
        sigmaP = 2 * Math.PI;
    while (Math.abs(sigma - sigmaP) > 1e-12) {
        var cos2SigmaM = Math.cos(2*sigma1 + sigma),
            sinSigma = Math.sin(sigma),
            cosSigma = Math.cos(sigma),
            deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
            sigmaP = sigma;
        sigma = s / (b * A) + deltaSigma;
    };
    var tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1,
        lat2 = Math.atan2(sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1, (1 - f) * Math.sqrt(sinAlpha * sinAlpha + tmp * tmp)),
        lambda = Math.atan2(sinSigma * sinAlpha1, cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1),
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha)),
        L = lambda - (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM))),
        revAz = Math.atan2(sinAlpha, -tmp); // final bearing
    return [toDeg(lat2), lon1 + toDeg(L)];
};

function gpxOneCircle(lat, lng, radius)
{
    var gpxResult = "";
    var pt = destVincenty(lat,lng,0,radius);
    gpxResult += '<rtept lat="' + pt[0] + '" lon="' + pt[1] + '" />';
    for (var d=0; d != 721; ++d)
    {
        pt = destVincenty(lat,lng,d*0.5,radius);
        gpxResult += '<rtept lat="' + pt[0] + '" lon="' + pt[1] + '" />';

    }
    return gpxResult;
}

function gpxCircles() {
    var gpxResult = "";
    gpxResult += "\n<rte><name>Cercle de 1km</name>";
    gpxResult += gpxOneCircle(window.center.lat, window.center.lng, 1000);
    gpxResult += "</rte>";
    gpxResult += "\n<rte><name>Cercle de 10km</name>";
    gpxResult += gpxOneCircle(window.center.lat, window.center.lng, 10000);
    gpxResult += "</rte>";
    return gpxResult;
}
function gpxForLoop() {
    var gpx = "";
    
    gpx += gpxHeader;
    if (window.ClickCircularPaths.length != 0) { 
       gpx += gpxBody;
    }
    gpx += gpxCircles();
    gpx += gpxFooter;

    return gpx;
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
   
    const square = walkCircle.getBounds();
    const bboxString = square.getSouth() + "," + square.getWest() + ", " + square.getNorth() + "," + square.getEast();

    const query = "[out:json][timeout:25];(way[\"highway\"][\"highway\"!~\"^(motorway|construction|trunk|trunk_link|motorway_link)\"][\"foot\"][\"foot\"!~\"^(yes|designated|permissive)\"](" + bboxString + ");\
    way[\"highway\"][\"highway\"!~\"^(motorway|construction|trunk|trunk_link|motorway_link)\"][\"access\"!~\"^(no|service|private)\"][!\"foot\"](" + bboxString + ");\
    way[\"pathway\"](" + bboxString + "););out;>;out skel qt;";
    
    const server = "https://overpass.kumi.systems/api/interpreter";

    const request = new XMLHttpRequest();
    
    
    request.open('GET', server + "?data=" + encodeURIComponent(query), true);
    request.setRequestHeader('Content-type', 'application/json');

    request.addEventListener("load", function() {

        const idProcess = window.idProcess;

        const response = JSON.parse(request.response);

        gpxBody = "";

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

                // draw the shortest path and generate GPX
                if (shortestPath.edges.length)
                {
                    gpxBody += "<rte><name>Accès à la boucle</name>";
                    gpxBody += '<rtept lat="' + shortestPath.nodes[shortestPath.edges[0][0]].lat + '" lon="' + shortestPath.nodes[shortestPath.edges[0][0]].lng + '" />';
                            for(var e = 0; e != shortestPath.edges.length; ++e) {
                    gpxBody += '<rtept lat="' + shortestPath.nodes[shortestPath.edges[e][1]].lat + '" lon="' + shortestPath.nodes[shortestPath.edges[e][1]].lng + '" />';
                    window.ClickCircularPaths.push(L.polyline(toPolyLine(shortestPath.edges[e], shortestPath.nodes), { color: '#2ab50a'}).addTo(map));
                            }
                    gpxBody += "</rte>\n";
                }

                var neighbours = {};

                // draw the contour and generate GPX
                for(var e = 0; e != graph.edges.length; ++e) {
                    window.ClickCircularPaths.push(L.polyline(toPolyLine(graph.edges[e], graph.nodes), { color: '#0060f0'}).addTo(map));

                    if (! neighbours[graph.edges[e][0]]) neighbours[graph.edges[e][0]] = [];
                    neighbours[graph.edges[e][0]].push([graph.edges[e][1]]);
                    if (! neighbours[graph.edges[e][1]]) neighbours[graph.edges[e][1]] = [];
                    neighbours[graph.edges[e][1]].push([graph.edges[e][0]]);
                }

                var neighboursVisited = {};
                var currentNode;
                var nextNode;

                if (shortestPath.edges[0])
                {
                    currentNode = shortestPath.edges[0][0];
                }
                else
                {
                    currentNode = graph.edges[0][0];
                }

                gpxBody += "<rte><name>Boucle de marche</name>";
                gpxBody += '<rtept lat="' + graph.nodes[currentNode].lat + '" lon="' + graph.nodes[currentNode].lng + '" />';

                nextNode = neighbours[currentNode].pop();
                while (nextNode)
                {
                    gpxBody += '<rtept lat="' + graph.nodes[nextNode].lat + '" lon="' + graph.nodes[nextNode].lng + '" />';
                    if (! neighboursVisited[currentNode]) neighboursVisited[currentNode] = [];
                    neighboursVisited[currentNode].push(nextNode);
                    if (! neighboursVisited[nextNode]) neighboursVisited[nextNode] = [];
                    neighboursVisited[nextNode].push(currentNode);
                    currentNode = nextNode;
                    nextNode = neighbours[currentNode].pop();
                    while (neighboursVisited[nextNode] && neighboursVisited[nextNode].includes(currentNode)) nextNode = neighbours[currentNode].pop();
                }
                gpxBody += "</rte>\n";

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
