function log(x) {
    console.debug(x);
}


HalfEdge = Object();

HalfEdge.buildHalfEdges = function(geometry) {
    // loop through faces and build a half-edge for each edge of each face,
    // containing one vertex, one face and some other stuff

    var i, j, A, B, count, f,
        faces = geometry.faces,
        edgeCount = faces.length * 3,
        keyEdge = Array(geometry.vertices.length),
        edgeVert = Array(edgeCount),
        edgeFace = Array(edgeCount),
        edgeNext = Array(edgeCount),
        edgePair = Array(edgeCount);

    // for each vertex, create a list of vertex pairs for each neighboring triangle
    for(i=0, j=0, count=faces.length; i<count; i++) {
        f = geometry.faces[i];
        if(f.a === f.b) { // no degenerates!
            j += 3;
            throw 'Encountered degenerate face ['+i+']: '+f.a +', '+f.b+', '+f.c+", ain't nobody got time for that!";
            // continue;
            // this is probably going to leave undefined entries in the arrays
            // and break some stuff looping over them
        }

          //         A    Edge shown has vert = B, face = ABC,
          //        / \     and might be the key edge of A
          //       /   \
          //      / /   \
          //     / /     \
          //    / /-      \
          //   /           \
          //  B ----------- C

        keyEdge[f.c] = j;   // the key edge is pointing to the next vertex
        edgeVert[j] = f.a;
        edgeFace[j] = i;
        edgeNext[j] = ++j;

        keyEdge[f.a] = j;
        edgeVert[j] = f.b;
        edgeFace[j] = i;
        edgeNext[j] = ++j;

        keyEdge[f.b] = j;
        edgeVert[j] = f.c;
        edgeFace[j] = i;
        edgeNext[j] = ++j - 3;
    }

    // now find the edge pairs
    outerLoop:
    for(i=0, count=edgePair.length; i<count; i++) {
        if(edgePair[i] !== undefined)
            continue;

        A = edgeVert[i]; // primary vertex of this edge
        B = edgeVert[edgeNext[edgeNext[i]]]; // opposite vertex

        for(j=0; j<count; j++) {

            if(edgeVert[j] === B && edgeVert[edgeNext[edgeNext[j]]] === A) {
                edgePair[i] = j;
                edgePair[j] = i;
                continue outerLoop;
            }
        }
        throw "Could not find opposite vertex! :(";
    }


    geometry.freeEdges = Array();
    geometry.freeFaces = Array();
    geometry.freeVerts = Array();

    geometry.edges = Object();

    geometry.edges.vert = edgeVert;
    geometry.edges.face = edgeFace;
    geometry.edges.next = edgeNext;
    geometry.edges.pair = edgePair;
    geometry.edges.key = keyEdge;
    geometry.edges.lengthSq = Array(edgeCount);
    geometry.edges.deleted = Array(edgeCount);
    for(i=0; i<edgeCount; i++)
        geometry.edges.deleted[i] = false;
};

HalfEdge.checkEdgePairs = function(geometry) {
    var i, count,
    edgePair = geometry.edges.pair,
    deleted = geometry.edges.deleted;

    for(i=0, count=edgePair.length; i<count; i++) {
        if(deleted[i])
            continue;
        if(edgePair[edgePair[i]] !== i)
            throw "whoopsadaisy";
    }
};

HalfEdge.findLongestCycle = function(geometry) {
    var keyEdge = geometry.edges.key,
        edgeNext = geometry.edges.next,
        edgePair = geometry.edges.pair,
        // lengthSq = geometry.edges.lengthSq,
        // deleted = geometry.edges.deleted,
        i=0, j, count = keyEdge.length,
        edge, firstEdge, largest = 0;

    // loop over vertices
    for(i=0; i<count; i++) {
        j = 0;
        // loop around spokes
        edge = firstEdge = keyEdge[i];
        while(true) {
            j++;
            edge = edgeNext[edgePair[edge]];

            if(edge === firstEdge)
                break;
        }

        if(j > largest)
            largest = j;
    }
    return largest;
}

HalfEdge.computeEdgeLengths = function(geometry) {
    var vertices = geometry.vertices,
        edgeVert = geometry.edges.vert,
        edgePair = geometry.edges.pair,
        lengthSq = geometry.edges.lengthSq,
        deleted = geometry.edges.deleted,
        i, A, B, count = lengthSq.length;

    // first clear the old lengths (except deleted edges)
    for(i=0; i<count; i++) {
        if(!deleted[i])
            lengthSq[i] = -1;
    }

    // now calculate new ones, skipping any that were
    // already calculated for their pair (or are deleted)
    for(i=0; i<count; i++) {
        if(lengthSq[i] === -1) {
            A = vertices[edgeVert[i]];
            B = vertices[edgeVert[edgePair[i]]];

            // set length for both half-edges
            // if(A.x === -1000 || B.x === -1000)
            //     lengthSq[i] = lengthSq[edgePair[i]] = null;
            // else
            lengthSq[edgePair[i]] = lengthSq[i] = A.distanceToSquared(B);
        }
    }
}

HalfEdge.computeVertEdgeLengths = function(geometry, vertIndex) {
    HalfEdge.checkEdgePairs(geometry);

    var edges = geometry.edges,
        firstEdge = edges.key[vertIndex],
        edge = firstEdge,
        vertices = geometry.vertices,
        X = vertices[vertIndex], Y, dist, n = 0;

    while(true) {
        Y = vertices[edges.vert[edge]];
        dist = X.distanceToSquared(Y);
        edges.lengthSq[edge] = dist;
        edge = edges.pair[edge];
        edges.lengthSq[edge] = dist;

        edge = edges.next[edge];

        if(edge === firstEdge)
            break;
        if(++n > 20)
            throw("OH GOD. THERE'S SO MUCH BLOOD. MAKE IT STOP!!!1");
    }
};

// this.centroid = new THREE.Vector3(0,0,0)
// this.volume = 0;

// nifty little algorithm that sums positive and negative volumes
// of tetrahedrons from the origin to each face.
// also calculates the center of mass of the mesh
HalfEdge.calculateVolume = function(geometry) {
    var faces = geometry.faces,
        vertices = geometry.vertices,
        face,
        volume = 0,
        // centroid,
        A, B, C;
    // this.centroid.set(0,0,0);

    for(var i=0, count=faces.length; i<count; i++) {
        face = faces[i];
        A = vertices[face.a];
        B = vertices[face.b];
        C = vertices[face.c];

        // skip degenerate faces
        if(A === B || A === C || B === C)
            continue;

        volume += A.dot(B.clone().cross(C));

        // centroid = A.clone().add(B).add(C).multiplyScalar(volume / 4);
        // this.centroid.add(centroid);
    }

    // this.centroid.divideScalar(this.volume);
    return volume;
};

/*// calculate the mean velocity of the mesh, and the mean solid-body rotation
MeshUtils.calculateSolidMotion = function() {
    var realVertices = 0, vertex;

    var v = new THREE.Vector3(0,0,0);
    var omega = new THREE.Vector3(0,0,0);

    for(i=0, count=vertices.length; i<count; i++) {
        if(vertices[i].x === -1000)
            continue;
        realVertices++; // count only vertices in use currently

        v.add(velocity[i]);

        omega.add(vertices[i].clone().sub(this.centroid).cross(velocity[i]));

    }

    v.divideScalar(realVertices);
    omega.divideScalar(realVertices);

}*/

HalfEdge.removeEdge = function(geometry, i) {
    var edges = geometry.edges;
    // if(edges.vert[i] === null) {
        // debugger;
    //     return;
    // }
    // log('removing edge '+i);



    // for(var j=0, count=edges.key.length; j<count; j++) {
    //     if(edges.key[j] == i)
    //         throw("uh-oh");
    // }

    // make sure we relace it if this edge was referenced by its vertex
    // if(edges.key[edges.vert[edges.pair[i]]] === i) {
    //     throw("hey you noob, this edge is still the key of vert " + edges.vert[edges.pair[i]]);
        // edges.key[edges.vert[i]] = edges.pair[edges.next[i]];
        // console.log('changing keyEdge['+this.edgeVert[i]+']  from '+i+' to '+this.edgePair[this.edgeNext[i]]);
    // }

    edges.deleted[i] = true;
    edges.vert[i] = null;
    edges.face[i] = null;
    edges.next[i] = null;
    edges.pair[i] = null;
    edges.lengthSq[i] = null;
    geometry.freeEdges.push(i);
};

HalfEdge.removeFace = function(geometry, i) {
    // log('removing face '+i);
    var face = geometry.faces[i];
    face.a = face.b = face.c = 0;
    geometry.freeFaces.push(i);
}

HalfEdge.removeVert = function(geometry, i) {
    // log('removing vert '+i);
    var vert = geometry.vertices[i];
    vert.x = vert.y = vert.z = 0;
    geometry.edges.key[i] = null;
    geometry.freeVerts.push(i);
}

HalfEdge.mergeEdge = function(geometry, AX) {

    // log('Merging edge '+AX)

    var i, count, edge, lastEdge, face,
        vertices = geometry.vertices,
        faces = geometry.faces,
        edgeVert = geometry.edges.vert,
        keyEdge = geometry.edges.key,
        edgePair = geometry.edges.pair,
        edgeNext = geometry.edges.next,
        edgeFace = geometry.edges.face;

    // if(this.edgeLength[AX] === null)
    //     throw('dangit');

    // -*-------L------L2---   we want to delete X and reconnect L2, R2, etc to A
    // / \     /`\     / \
    //    \   /```\   /   \    move A half way towards X and average their velocities
    //     \ /`````\ /     \
    // -----A-------X-------*
    //     / \`````/ \     /
    //    /   \```/   \   /
    // \ /     \`/     \ /
    // -*-------R-------R2---

    var XA = edgePair[AX],
        XL = edgeNext[AX], LA = edgeNext[XL],
        AR = edgeNext[XA], RX = edgeNext[AR],
        XR = edgePair[RX], LX = edgePair[XL],
        RR2 = edgeNext[XR], R2X = edgeNext[RR2],
        XL2 = edgeNext[LX], L2L = edgeNext[XL2];

    var A = edgeVert[XA],
        R = edgeVert[XR],
        L = edgeVert[XL],
        L2 = edgeVert[XL2],
        R2 = edgeVert[RR2],
        X = edgeVert[AX];

    // make sure vertices A, R and L don't end up without key edges
    keyEdge[A] = AR;
    keyEdge[L] = LA;
    keyEdge[R] = RR2;

    // replace faces that are about to be deleted
    edgeFace[LA] = edgeFace[LX];
    edgeFace[AR] = edgeFace[XR];

    // loop through the "spokes" of vertex X
    var n = 2, edge = LX;
    while(edge !== RX) {
        // replace vertex X with A in the face
        face = faces[edgeFace[edge]];

        if(face.a === X)
            face.a = A;
        else if(face.b === X)
            face.b = A;
        else if(face.c === X)
            face.c = A;
        else
            throw 'uh-oh: '+face.a+', '+face.b+', '+face.c;

        // update the vertex of each inward edge
        edgeVert[edge] = A;

        edge = edgePair[edgeNext[edge]];

        if(n++ > 20)
            throw 'whoops, infinite loop on aisle 3!';
    }

    // repair the next-edge relationships

    // ---------------L--     degenerate case where X
    // \          _-'/|\      has only 3 neighbors
    //  \      _-'  / | \
    //   \  _-'    /  |
    // ---A-------X   |
    //   / `-_     \  |
    //  /     `-.   \ | /
    // /         `-._\|/
    // ---------------R---
    // if(R === L2) { // also R2 === L
    if(RR2 === L2L) {
        edgeNext[LA] = AR;
        edgeNext[L2L] = LA;

    } else {
        edgeNext[LA] = XL2;
        edgeNext[R2X] = AR;
        edgeNext[L2L] = LA;
    }
    edgeNext[AR] = RR2;

    // delete stuff
    this.removeFace(geometry, edgeFace[AX]);
    this.removeFace(geometry, edgeFace[XA]);

    this.removeEdge(geometry, AX);
    this.removeEdge(geometry, XA);
    this.removeEdge(geometry, LX);
    this.removeEdge(geometry, XL);
    this.removeEdge(geometry, RX);
    this.removeEdge(geometry, XR);


    this.removeVert(geometry, X);

    /*HalfEdge.checkEdgePairs(geometry)

    // check if anything didn't get deleted good
    var removedEdges = [AX, XA, LX, XL, RX, XR];
    var deleted = geometry.edges.deleted;


    for(i=0, count=keyEdge.length; i<count; i++) {
        if(i === X || keyEdge[i] === null)
            continue;
        if(edgeVert[edgePair[keyEdge[i]]] !== i)
            throw("uh-oh");
    }


    for(i=0, count=edgeVert.length; i<count; i++) {
        if(deleted[i])
            continue;

        if(deleted[edgeNext[i]])
            throw("UGHHHH");
        if(deleted[edgePair[i]])
            throw("UGHHHH");
    }

    for(i=0, count=keyEdge.length; i<count; i++) {
        if(geometry.freeVerts.indexOf(i) > -1)
            continue;
        if(deleted[keyEdge[i]])
            throw("UGHHHH");
    }


    for(i=0, count=edgeVert.length; i<count; i++) {
        if(edgeVert[i] == X)
            throw("UGHHHH");
    }*/
};

HalfEdge.splitEdge = function(AB) {
    var freeFaces = this.freeFaces,
        freeVerts = this.freeVerts,
        freeEdges = this.freeEdges,
        edgeVert = this.edgeVert,
        edgePair = this.edgePair,
        edgeNext = this.edgeNext,
        edgeFace = this.edgeFace,
        edgeLength = this.edgeLength,
        vertices = this.geometry.vertices,
        face,
        faces = this.geometry.faces;

    if(freeFaces.length < 2 || freeVerts.length === 0 || freeEdges.length < 6)
        return null; // this edge lives to see another day....for now

    log('splitting edge '+AB+' ['+edgeVert[AB]+'->'+edgeVert[edgePair[AB]]+']');

    //           B
    //          /|\
    //         / | \
    //        /  |  \
    //       /   |   \
    //      L----X----R
    //       \   |   /
    //        \  |  /
    //         \ | /
    //          \|/
    //           A

    var BA = edgePair[AB],
        AR = edgeNext[BA],
        RB = edgeNext[AR],
        BL = edgeNext[AB],
        LA = edgeNext[BL],

        A = edgeVert[BA],
        B = edgeVert[AB],
        R = edgeVert[AR],
        L = edgeVert[BL],
        X = freeVerts.pop(),

        AXL = edgeFace[AB], // reusing ABL
        XBL = edgeFace[BA], // reusing BAR
        ARX = freeFaces.pop(),
        XRB = freeFaces.pop();

    this.removeEdge(AB);
    this.removeEdge(BA);

    var AX = freeEdges.pop(),
        XA = freeEdges.pop(),
        XB = freeEdges.pop(),
        BX = freeEdges.pop(),
        XL = freeEdges.pop(),
        LX = freeEdges.pop(),
        XR = freeEdges.pop(),
        RX = freeEdges.pop();

    face = faces[AXL], face.a = A, face.b = X, face.c = L,
    face = faces[XBL], face.a = X, face.b = B, face.c = L,
    face = faces[ARX], face.a = A, face.b = R, face.c = X,
    face = faces[XRB], face.a = X, face.b = R, face.c = B,

    edgeNext[AX] = XL, edgeVert[AX] = X, edgeFace[AX] = AXL, edgePair[AX] = XA,
    edgeNext[XA] = AR, edgeVert[XA] = A, edgeFace[XA] = ARX, edgePair[XA] = AX,
    edgeNext[XB] = BL, edgeVert[XB] = B, edgeFace[XB] = XBL, edgePair[XB] = BX,
    edgeNext[BX] = XR, edgeVert[BX] = X, edgeFace[BX] = XRB, edgePair[BX] = XB,
    edgeNext[XL] = LA, edgeVert[XL] = L, edgeFace[XL] = AXL, edgePair[XL] = LX,
    edgeNext[LX] = XB, edgeVert[LX] = X, edgeFace[LX] = XBL, edgePair[LX] = XL,
    edgeNext[XR] = RB, edgeVert[XR] = R, edgeFace[XR] = XRB, edgePair[XR] = RX,
    edgeNext[RX] = XA, edgeVert[RX] = X, edgeFace[RX] = ARX, edgePair[RX] = XR;

    // reset key edges
    this.keyEdge[X] = AX;
    this.keyEdge[A] = XA;
    this.keyEdge[B] = XB;

    return X;
};
