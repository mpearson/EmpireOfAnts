
// THREE.Vector3.prototype.edge = null;
// THREE.Face3.prototype.edge = null;



function log(x) {
    // console.debug(x);
}

// water drop class
function WaterDrop() {
    // frames between mesh optimizations
    this.cleanupInterval = 4;

    // number of evolutions per frame
    this.iterations = 1;

    //
    this.timeFactor = 1;

    var size = 1,
        aspectRatio = 1,
        detail = 16,
        segmentSize = size/detail,
        noiseFactor = 0.0 * segmentSize;

    // segment length threshholds as a factor of the "target" length
    this.mergeThreshold = Math.pow(0.8 * segmentSize, 2);
    this.splitThreshold = Math.pow(1.8 * segmentSize, 2);

    // actual surface tension strength
    this.forceFactor = (0.04 / segmentSize) * this.timeFactor / this.iterations;

    // decays velocity with every frame
    this.damping = -1.0 * segmentSize * this.timeFactor / this.iterations;
    // this.normalDamping = -2.5 * segmentSize * this.timeFactor / this.iterations;
    // this.tangentDamping = -5.0 * segmentSize * this.timeFactor / this.iterations;

    // local velocity "blurring"...bah, doesnt' work
    this.velocityBleed = 0.0;

    // local volume correction: factor used to apply negative 'reaction'
    // acceleration to surrounding vertices
    this.lvcFactor = -1.0;

    this.geometry = new THREE.BoxGeometry(
        size * aspectRatio,
        size,
        size,
        detail * aspectRatio,
        detail,
        detail);
    // this.geometry = new THREE.TorusGeometry(0.8, 0.2, 16, 50);
    // this.geometry = new THREE.IcosahedronGeometry(1,3);
    this.geometry.mergeVertices();

    var vertices = this.geometry.vertices;
    var faces = this.geometry.faces;
    var vertexCount = vertices.length;
    var i, count;

    this.material = new THREE.MeshBasicMaterial({color: 0x2a87ff, wireframe: true });
    // this.material = new THREE.MeshBasicMaterial({color: 0x4499ff, wireframe: true });

    // this.material = new THREE.MeshBasicMaterial({color: 0x4499ff});
    // this.material = new THREE.MeshBasicMaterial({color: 0x447733, wireframe: true });
    // this.material = new THREE.MeshBasicMaterial({color: 0x87c9ff, wireframe: true });
    // this.material = new THREE.MeshBasicMaterial({color: 0xff00ff, wireframe: true });
    this.mesh = new THREE.Mesh(this.geometry, this.material);

    this.addNoise = function() {
        // introduce some noise
        for(i=0, count=vertices.length; i<count; i++) {

            // this.geometry.vertices[i].x *= 0.5;

            vertices[i].y += noiseFactor * (Math.random() - 0.5);
            vertices[i].x += noiseFactor * (Math.random() - 0.5);
            vertices[i].z += noiseFactor * (Math.random() - 0.5);
        }
    };

    this.addNoise();

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
        var A, B, C, totalLength = 0, f, AB, count = 1;

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
            C = vertices[edgeVert[edgeNext[edge]]];
            // AB = B.clone().sub(A);
            f = C.clone().sub(B).multiplyScalar(0.5).add(B).sub(A);
            // totalLength += AB.length();

            netForce.add(f);

            edge = edgeNext[edgePair[edge]];
            if(edge === firstEdge)
                break;
            count++;
        }

        netForce.multiplyScalar(this.forceFactor / count);
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
        if(this.velocityBleed == 0)
            return;
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
        v.multiplyScalar(1 + this.damping);

        // if(this.vertexNormals[i] !== undefined) {
        //     var vTangent = v.clone().projectOnPlane(this.vertexNormals[i]);
        //     v.sub(vTangent);    // v is now just the normal component

        //     v.multiplyScalar(this.normalDamping);
        //     v.add(vTangent.multiplyScalar(this.tangentDamping));
        // }
    }

    this.startingVolume = HalfEdge.calculateVolume(this.geometry);

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
        var iteration, i, count, AB, A, B, X, vA, vB, len;
        var vertices = this.geometry.vertices,
            velocity = this.velocity,
            keyEdge = this.geometry.edges.key,
            deleted = this.geometry.edges.deleted,
            edgeVert = this.geometry.edges.vert,
            edgePair = this.geometry.edges.pair,
            edgeLength = this.geometry.edges.lengthSq;



        for(iteration=0; iteration<this.iterations; iteration++) {
            for(i=0, count=vertices.length; i<count; i++) {
                if(keyEdge[i] === null)
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

                for(AB=0, count=edgeLength.length; AB<count; AB++) {
                    if(deleted[AB])
                        continue;

                    len = edgeLength[AB];

                    A = edgeVert[edgePair[AB]];
                    B = edgeVert[AB];

                    if(len < this.mergeThreshold) {
                        if(merged.indexOf(A) > -1) {// || merged.indexOf(B) > -1) {
                            // log('skipping edge '+AB);
                            continue;
                        }

                        vA = velocity[A];
                        vB = velocity[B];

                        // average the velocities
                        vA.add(vB).multiplyScalar(0.5);
                        // vB.set(0, 0, 0);

                        // move B half way to A
                        vertices[A].add(vertices[B]).multiplyScalar(0.5);

                        HalfEdge.mergeEdge(this.geometry, AB);
                        // HalfEdge.computeEdgeLengths(this.geometry);
                        HalfEdge.computeVertEdgeLengths(this.geometry, A);

                        merged.push(A);

                    } else if(len > this.splitThreshold) {
                        if(merged.indexOf(A) > -1 || merged.indexOf(B) > -1) {
                            // log('skipping edge '+AB);
                            continue;
                        }
                        X = HalfEdge.splitEdge(this.geometry, AB);
                        if(X !== null) {
                            // X has average the velocities of A and B
                            velocity[X] = velocity[B].clone().add(velocity[A]).multiplyScalar(0.5);

                            // set X to half way half way between A and B
                            vertices[X] = vertices[B].clone().add(vertices[A]).multiplyScalar(0.5);

                            HalfEdge.computeVertEdgeLengths(this.geometry, X);
                            merged.push(A);
                            merged.push(B);
                            merged.push(X);
                        }
                    }
                }
            }
        }
        this.geometry.computeFaceNormals();
        this.geometry.computeVertexNormals();
    };

    this.frameCount = 0;
}
