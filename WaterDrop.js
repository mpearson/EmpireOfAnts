
// THREE.Vector3.prototype.edge = null;
// THREE.Face3.prototype.edge = null;



function log(x) {
    console.debug(x);
}

// water drop class
function WaterDrop() {
    this.cleanupInterval = 4; // frames between mesh optimizations

    this.size = 1;
    this.detail = 10;
    this.mergeThreshold = Math.pow(0.9 * this.size / this.detail, 2);
    this.splitThreshold = Math.pow(1.9 * this.size / this.detail, 2);
    this.aspectRatio = 2;

    this.forceFactor = 0.01;
    // this.normalDamping = 0.99;
    // this.tangentDamping = 0.9;
    this.damping = 0.99;
    this.velocityBleed = 0;//.05;
    this.iterations = 1;
    this.noiseFactor = 0.0;
    this.lvcFactor = -1;

    this.geometry = new THREE.BoxGeometry(
        this.size*this.aspectRatio,
        this.size,
        this.size,
        this.detail*this.aspectRatio,
        this.detail,
        this.detail);
    // this.geometry = new THREE.TorusGeometry(0.8, 0.2, 16, 50);
    // this.geometry = new THREE.IcosahedronGeometry(1,3);
    this.geometry.mergeVertices();


    var vertices = this.geometry.vertices;
    var faces = this.geometry.faces;
    var vertexCount = vertices.length;
    var i, count;

    this.material = new THREE.MeshBasicMaterial({color: 0x4499ff, wireframe: true });
    // this.material = new THREE.MeshBasicMaterial({color: 0x447733, wireframe: true });
    // this.material = new THREE.MeshBasicMaterial({color: 0x87c9ff, wireframe: true });
    // this.material = new THREE.MeshBasicMaterial({color: 0xff00ff, wireframe: true });
    this.mesh = new THREE.Mesh(this.geometry, this.material);

    this.addNoise = function(magnitude) {
        // introduce some noise
        for(i=0, count=vertices.length; i<count; i++) {

            // this.geometry.vertices[i].x *= 0.5;

            vertices[i].y += magnitude * (Math.random() - 0.5);
            vertices[i].x += magnitude * (Math.random() - 0.5);
            vertices[i].z += magnitude * (Math.random() - 0.5);
        }
    };

    this.neighborMap = null;
    this.vertexNormals = null;
    this.faceMap = null;

    this.velocity = new Array(vertices.length);
    var dv = this.velocity;

    for(i=0, count=vertices.length; i<count; i++)
        dv[i] = new THREE.Vector3(0,0,0);



    this.mcfMaterial = new THREE.LineBasicMaterial({color: 0x00ff00});
    this.mcfHelperGeom = new THREE.Geometry();
    this.mcfHelperObject = new THREE.Line(this.mcfHelperGeom, this.mcfMaterial);
    this.mcfHelperObject.type = THREE.LinePieces;

    this.velocityMaterial = new THREE.LineBasicMaterial({color: 0xff0000});
    this.velocityHelperGeom = new THREE.Geometry();
    this.velocityHelperObject = new THREE.Line(this.velocityHelperGeom, this.velocityMaterial);
    this.velocityHelperObject.type = THREE.LinePieces;

    this.mcfHelperGeom.vertices = new Array(vertices.length * 2);
    this.velocityHelperGeom.vertices = new Array(vertices.length * 2);
    for(i=0, count=vertices.length; i<count; i++) {
        this.mcfHelperGeom.vertices[i] = new THREE.Vector3(0,0,0);
        this.velocityHelperGeom.vertices[i] = new THREE.Vector3(0,0,0);
    }

    // calculate the weighted sum of forces from each neighboring vertex
    // ...the "surface tension" if you will
    this.applySurfaceTension = function(i) {
        var vertices = this.geometry.vertices,
            velocity = this.velocity,
            keyEdge = this.geometry.edges.key,
            edgeVert = this.geometry.edges.vert,
            edgeNext = this.geometry.edges.next,
            edgePair = this.geometry.edges.pair;

        var netForce = new THREE.Vector3(0,0,0);
        var A, B, BA, count = 1;

        A = vertices[i];

        // ------B-------B-----
        //      / \     / \
        // \   /   \   /   \
        //  \ /     \ /     \ /
        // --B-------A-------B
        //  / \     / \     / \
        // /   \   /   \   /
        //      \ /     \ /
        // ------B-------B----

        // loop around the neighborhood! weeee!
        var firstEdge = keyEdge[i];
        var edge = firstEdge;
        var iters = 0;
        while(true) {
            if(iters++ > 20)
                throw('dangit');
            B = vertices[edgeVert[edge]];
            netForce.add(B).sub(A);

            edge = edgeNext[edgePair[edge]];
            if(edge === firstEdge)
                break;
            count++;
        }

        netForce.multiplyScalar(this.forceFactor);
        velocity[i].add(netForce);

        netForce.multiplyScalar(this.lvcFactor / count);

        // around we go again!
        edge = firstEdge;
        while(true) {
            velocity[edgeVert[edge]].add(netForce);

            edge = edgeNext[edgePair[edge]];
            if(edge === firstEdge)
                break;
        }
    };

    this.localAverage = function(i) {
        var firstEdge = this.geometry.edges.key[i],
            velocity = this.velocity,
            edgeVert = this.geometry.edges.vert,
            edgeNext = this.geometry.edges.next,
            edgePair = this.geometry.edges.pair;
            edge = firstEdge,
            count = 0,
            average = new THREE.Vector3(0,0,0);

        while(true) {
            count++;
            average.add(velocity[edgeVert[edge]]);

            edge = edgeNext[edgePair[edge]];
            if(edge === firstEdge)
                break;
        }
        average.multiplyScalar(this.velocityBleed/count);

        velocity[i].multiplyScalar(1-(this.velocityBleed)).add(average);
    };


    this.applyVelocity = function(i) {
        // apply dV to vertex
        this.localAverage(i);
        var v = this.velocity[i];

        this.geometry.vertices[i].add(v);
        v.multiplyScalar(this.damping);

        // if(this.vertexNormals[i] !== undefined) {
        //  var vTangent = v.clone().projectOnPlane(this.vertexNormals[i]);
        //  v.sub(vTangent);    // v is now just the normal component

        //  v.multiplyScalar(this.normalDamping);
        //  v.add(vTangent.multiplyScalar(this.tangentDamping));
        // }
    }

    this.setStartingVolume = function() {
        this.startingVolume = HalfEdge.calculateVolume(this.geometry);
    }

    this.setStartingVolume();

    this.correctGlobalVolume = function() {
        this.volume = HalfEdge.calculateVolume(this.geometry);
        // this.calculateVolume();

        var correction = Math.pow(this.startingVolume/this.volume, 1/3),
            vert,
            deleted = this.geometry.edges.deleted,
            vertices = this.geometry.vertices;

        for(var i=0, count = vertices.length; i<count; i++) {
            if(deleted[i])    // skip deleted vertices
                continue;
            vert = vertices[i];
            vert.x *= correction;
            vert.y *= correction;
            vert.z *= correction;
        }
    };


    this.evolve = function() {
        var i, count, A, B, X, vA, vB, len;
        var vertices = this.geometry.vertices,
            velocity = this.velocity,
            vertEdge = this.geometry.edges.key,
            deleted = this.geometry.edges.deleted,
            edgeVert = this.geometry.edges.vert,
            edgePair = this.geometry.edges.pair,
            edgeLength = this.geometry.edges.lengthSq;

        for(i=0, count=vertices.length; i<count; i++) {
            if(vertEdge[i] === null)
                continue;

            this.applySurfaceTension(i);
            this.applyVelocity(i);
        }

        this.correctGlobalVolume();
        HalfEdge.computeEdgeLengths(this.geometry);

        if(this.frameCount++ === this.cleanupInterval) {
            this.frameCount = 0;
            log('>>>>>>>>>>>>>>>>mesh cleanup frame');

            var merged = Array();

            for(i=0, count=edgeLength.length; i<count; i++) {
                if(deleted[i])
                    continue;

                len = edgeLength[i];

                A = edgeVert[edgePair[i]];
                B = edgeVert[i];

                if(len < this.mergeThreshold) {
                    if(merged.indexOf(A) > -1 || merged.indexOf(B) > -1) {
                        // log('skipping edge '+i);
                        // continue;
                    }

                    vA = velocity[A];
                    vB = velocity[B];

                    // average the velocities
                    vA.add(vB).multiplyScalar(0.5);
                    // vB.set(0, 0, 0);

                    // move B half way to A
                    vertices[A].add(vertices[B]).multiplyScalar(0.5);

                    HalfEdge.mergeEdge(this.geometry, i);
                    // this.computeEdgeLengths();
                    HalfEdge.computeVertEdgeLengths(this.geometry, A);

                    merged.push(A);
                    merged.push(B);


                } else if(false && len > this.splitThreshold) {
                    if(merged.indexOf(A) > -1 || merged.indexOf(B) > -1) {
                        // log('skipping edge '+i);
                        continue;
                    }
                    X = this.splitEdge(i);
                    if(X !== null) {
                        // X has average the velocities of A and B
                        velocity[X] = velocity[B].clone().add(velocity[A]).multiplyScalar(0.5);

                        // set X to half way half way between A and B
                        vertices[X] = vertices[B].clone().add(vertices[A]).multiplyScalar(0.5);

                        this.computeVertexEdgeLengths(X);
                        merged.push(A);
                        merged.push(B);
                        merged.push(X);
                    }
                }
            }
        }

        this.geometry.computeFaceNormals();
        this.geometry.computeVertexNormals();
    };



    this.frameCount = 0;
}
