
// THREE.Vector3.prototype.edge = null;
// THREE.Face3.prototype.edge = null;


// water drop class
function WaterDrop() {


	this.size = 1;
	this.detail = 16;
	this.mergeThreshold = 0.02 * this.size / this.detail;
	this.splitThreshold = 0.06 * this.size / this.detail;
	this.aspectRatio = 1.5;


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

	this.material = new THREE.MeshBasicMaterial({color: 0x447733, wireframe: true });
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
	}

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

	this.sortNeighbors = function(pairs) {
		if(pairs.length < 2)
			return [];

		var i, j, count = pairs.length/2, count2 = pairs.length;
		var neighbors = Array(count);

		neighbors[0] = pairs[0];
		neighbors[1] = pairs[1];

		var position = 1;

		// search for a pair containing the current vertex, and add its comrade to the end of the array
		for(i=0; i < count; i++) {
			// we can skip the first 2 since they've already been added
			for(j = 2; j < count2; j += 2) {
				if(pairs[j] === neighbors[position])
					neighbors[++position] = pairs[j +1 ];
				else if(pairs[j + 1] === neighbors[position])
					neighbors[++position] = pairs[j];
				else
					continue;
				pairs[j] = pairs[j + 1] = null;	// once a vertex is added, we don't want to find it again
				break;
			}
		}
		return neighbors;
	};

	// build a data structure mapping a vertex to an ordered array of the surrounding vertices
	this.buildNeighborMap = function() {
		var vertices = this.geometry.vertices;
		var faces = this.geometry.faces;
		var i, j, count, count2, face, neighbors;

		var neighborMap = Array(count);
		var vertexNormals = Array(count);

		for(i=0, count=vertices.length; i<count; i++)
			neighborMap[i] = new Array();

		// for each vertex, create a list of vertex pairs for each neighboring triangle
		for(i=0, count=faces.length; i<count; i++) {
			face = faces[i];
			if(face.a === face.b) {	// no degenerates!
				// console.log('skipping '+face.a +','+face.b+','+face.c);
				continue;
			}

			vertexNormals[face.a] = face.vertexNormals[0];
			vertexNormals[face.b] = face.vertexNormals[1];
			vertexNormals[face.c] = face.vertexNormals[2];

			neighbors = neighborMap[face.a];
			if(neighbors.indexOf(face.b) < 0)
				neighbors.push(face.b);
			if(neighbors.indexOf(face.c) < 0)
				neighbors.push(face.c);
			neighbors = neighborMap[face.b];
			if(neighbors.indexOf(face.a) < 0)
				neighbors.push(face.a);
			if(neighbors.indexOf(face.c) < 0)
				neighbors.push(face.c);
			neighbors = neighborMap[face.c];
			if(neighbors.indexOf(face.b) < 0)
				neighbors.push(face.b);
			if(neighbors.indexOf(face.a) < 0)
				neighbors.push(face.a);
		}
		this.neighborMap = neighborMap;
		this.vertexNormals = vertexNormals;
	};

	this.buildHalfEdges = function() {
		// loop through faces and create a half-edge for each face edge,
		// containing one vertex and some other stuff

		var vertices = this.geometry.vertices,
			faces = this.geometry.faces,
			i, j, count, face;


		var edgeVert = Array(faces.length * 3),
			edgeFace = Array(faces.length * 3),
			edgeNext = Array(faces.length * 3),
			edgePair = Array(faces.length * 3),

			vertEdge = Array(vertices.length);


		var vertexEdges = Array(vertices.length);

		// for each vertex, create a list of vertex pairs for each neighboring triangle
		for(i=0, j=0, count=faces.length; i<count; i++) {
			face = faces[i];
			if(face.a === face.b) {	// no degenerates!
				j += 3;
				// console.log('skipping '+face.a +','+face.b+','+face.c);
				continue;
			}

			// vertEdge[face.a] = j;
			edgeVert[j] = face.a;
			edgeFace[j] = i;
			edgeNext[j] = j++;

			// vertEdge[face.b] = j;
			edgeVert[j] = face.b;
			edgeFace[j] = i;
			edgeNext[j] = j++;

			// vertEdge[face.c] = j;
			edgeVert[j] = face.c;
			edgeFace[j] = i;
			edgeNext[j] = j++ - 2;
		}

		var A, B;
		loop1:
		for(i=0, count=edgeVert.length; i<count; i++) {
			if(edgePair[i] !== undefined)
				continue;

			A = edgeVert[i]; // primary vertex of this edge
			B = edgeVert[edgeNext[edgeNext[j]]]; // opposite vertex

			for(j=0; j<count; j++) {
				A = edgeVert[j];
				B = edgeVert[edgeNext[edgeNext[j]]];

				if(edgeVert[j] === B && edgeVert[edgeNext[edgeNext[j]]] === A) {
					edgePair[i] = j;
					edgePair[j] = i;
					continue loop1;
				}
			}
			throw "could not find opposite vertex D:";
		}

		this.vertEdge = vertEdge;
		this.edgeVert = edgeVert;
		this.edgeFace = edgeFace;
		this.edgeNext = edgeNext;
		this.edgePair = edgePair;

	};

	this.findLongestEdge = function() {
		var i, count, A, B, dist,
			vertices = this.geometry.vertices,
			edgeVert = this.edgeVert,
			edgeNext = this.edgeNext,
			edgePair = this.edgePair,
			maxDist = 0;

		for(i=0, count=edgeVert.length; i<count; i++) {
			A = edgeVert[i];
			B = edgeVert[edgeNext[edgeNext[i]]];

			dist = vertices[A].distanceTo(vertices[B]);
			if(dist > maxDist)
				maxDist = dist;


		}
		console.log('max dist: '+maxDist);

	}



	this.centroid = new THREE.Vector3(0,0,0)
	this.volume = 0;

	// nifty little algorithm that sums positive and negative volumes
	// of tetrahedrons from the origin to each face.
	// also calculates the center of mass of the mesh
	this.calculateVolume = function() {
		var face, volume, centroid, A, B, C;

		this.volume = 0;
		this.centroid.set(0,0,0);

		for(i=0, count=faces.length; i<count; i++) {
			face = faces[i];
			A = vertices[face.a];
			B = vertices[face.b];
			C = vertices[face.c];

			// skip degenerate faces
			if(A === B || A === C || B === C)
				continue;

			volume = A.dot(B.clone().cross(C));
			this.volume += volume;

			centroid = A.clone().add(B).add(C).multiplyScalar(volume / 4);
			this.centroid.add(centroid);
		}

		this.centroid.divideScalar(this.volume);
	};

	// calculate the mean velocity of the mesh, and the mean solid-body rotation
	this.calculateSolidMotion = function() {
		var realVertices = 0, vertex;

		var v = new THREE.Vector3(0,0,0);
		var omega = new THREE.Vector3(0,0,0);

		for(i=0, count=vertices.length; i<count; i++) {
			if(vertices[i].x === -1000)
				continue;
			realVertices++;	// count only vertices in use currently

			v.add(velocity[i]);

			omega.add(vertices[i].clone().sub(this.centroid).cross(velocity[i]));

		}

		v.divideScalar(realVertices);
		omega.divideScalar(realVertices);

	}


	this.calculateVolume();
	this.startingVolume = this.volume;

	this.forceFactor = 0.08;
	this.normalDamping = 0.99;
	this.tangentDamping = 0.9;
	this.damping = 0.98;
	this.iterations = 1;
	this.noiseFactor = 0.0;
	this.lvcFactor = -0.99;

	// calculate the weighted sum of forces from each neighboring vertex
	// ...the "surface tension" if you will
	this.applySurfaceTension = function(i) {
		var vertices = this.geometry.vertices;
		var netForce = new THREE.Vector3(0,0,0);
		var neighbors, n, count, A = vertices[i], BA, f, totalLength = 0;

		neighbors = this.neighborMap[i];

		for(n=0, count=neighbors.length; n<count; n++) {
			BA = A.clone().sub(vertices[neighbors[n]]);
			// totalLength += BA.length();

			netForce.sub(BA);
		}
		netForce.multiplyScalar(this.forceFactor);
		this.velocity[i].add(netForce);

		netForce.multiplyScalar(this.lvcFactor / count);

		for(n=0, count=neighbors.length; n<count; n++)
			this.velocity[neighbors[n]].add(netForce);
			// this.velocity[neighbors[n]].add(this.vertexNormals[neighbors[n]]).setLength(volumeDelta.length());
	};


	this.applyVelocity = function(i) {
		// apply dV to vertex
		var v = this.velocity[i];

		this.geometry.vertices[i].add(v);
		v.multiplyScalar(this.damping);

		// if(this.vertexNormals[i] !== undefined) {
		// 	var vTangent = v.clone().projectOnPlane(this.vertexNormals[i]);
		// 	v.sub(vTangent);	// v is now just the normal component

		// 	v.multiplyScalar(this.normalDamping);
		// 	v.add(vTangent.multiplyScalar(this.tangentDamping));
		// }
	}

	// this.mcfHelperEnabled = true;




	this.evolve2 = function() {




	};



	this.frameCount = 0;


	this.evolve = function() {
		var vertices = this.geometry.vertices;
		var i, j, k, neighborCount, vertexCount = vertices.length;
		var A, B, v, f, f2, neighbors, neighborIndex;



		for(i=0; i<vertexCount; i++) {
			A = vertices[i];
			if(A.x == -1000)
				continue;

			// v = this.velocity[i];


			// if(this.mcfHelperEnabled) {
			// 	this.mcfHelperGeom.vertices[2 * i] = A.clone().multiplyScalar(0.01);
			// 	this.mcfHelperGeom.vertices[2 * i + 1] = v.clone().multiplyScalar(0.01);	// store initial velocity for later

				// this.velocityHelperGeom.vertices[2 * i] = A.clone();
				// this.velocityHelperGeom.vertices[2 * i + 1] = A.clone().add(v.clone().multiplyScalar(200));
			// }

			this.applySurfaceTension(i);

			this.applyVelocity(i);

		}

		var mergeTargets = [];
		var merged = [];
		var neighborCount;
		var BA;
		var vNormal;


		// console.debug("evolving!");
		// return;
		for(i=0; i<vertexCount; i++) {
			A = vertices[i];
			v = this.velocity[i];


			if(this.mcfHelperEnabled) {
				this.mcfHelperGeom.vertices[2 * i + 1] = this.mcfHelperGeom.vertices[2 * i + 1].sub(v);
			}

			// find the 2 faces faces sharing A and B (F1 and F2)
			// find the far vertices of F1 and F2 (C1, C2)
			// remove F1 and F2 from the faceMap of A, B, C1 and C2
			//

			// replace A with B in A's neighborMap and the faces in A's faceMap
			if(merged.indexOf(i) > -1)
				continue;

			neighbors = this.neighborMap[i];

			if(this.frameCount++ === 5) {
				this.frameCount = 0;

				// loop through the neighboring vertices to see if any of them are too close
				for(j=0, neighborCount=neighbors.length; j<neighborCount; j++) {
					neighborIndex = neighbors[j];



					// skip ones that have already been merged with another point
					if(merged.indexOf(neighborIndex) > -1)
						continue;

					B = vertices[neighborIndex];

					BA = A.clone().sub(B);

					var len = BA.lengthSq()

					if(len < this.mergeThreshold) {
						// we must now eliminate vertex A, and replace all references to it with B
						mergeTargets[i] = neighborIndex;
						merged.push(i);
						merged.push(neighborIndex);

						// average the velocities
						this.velocity[neighborIndex].add(v).multiplyScalar(0.5);
						v.a = v.b = v.c = 0;

						// move B half way to A
						B.add(BA.multiplyScalar(0.5));

						// this.replaceNeighbor(i, B);


						// break;
					} else if(len > this.splitThreshold) {



					}
				}
			}

			// if(this.mcfHelperEnabled) {
				// this.mcfHelperGeom.vertices[2 * i + 1] = this.mcfHelperGeom.vertices[2 * i + 1].multiplyScalar(-0.2).add(A);
			// }
		}

		// if(merged.length == 0)
		// 	return;

		var faces = this.geometry.faces;
		var vIndex, face, faceCount = faces.length;

		// search for vertex A in each face and replace it with B
		// starting at the end so we can delete things
		for(i=faceCount-1; i>=0; i--) {
			face = faces[i];

			if((vIndex = mergeTargets[face.a]) !== undefined)
				face.a = vIndex;
			if((vIndex = mergeTargets[face.b]) !== undefined)
				face.b = vIndex;
			if((vIndex = mergeTargets[face.c]) !== undefined)
				face.c = vIndex;

			// if this face is now degenerate, then erase it from this world
			if(face.a === face.b || face.a === face.c || face.b == face.c) {
				// console.log('removing degenerate face ' + i + '('+face.a+','+face.b+','+face.c+')');
				face.a = face.b = face.c = 0;
			}
		}

		this.calculateVolume();

		var correctionFactor = Math.pow(this.startingVolume/this.volume, 1/3);
		// if(this.frameCount == 2) {
		// 	console.log(volume);
		// 	console.log(correctionFactor);
		// }

		for(i=0; i<vertexCount; i++) {
			if(mergeTargets[i] === undefined) {
				vertices[i].x *= correctionFactor;
				vertices[i].y *= correctionFactor;
				vertices[i].z *= correctionFactor;
			} else {
				// console.log('removing vertex ' + i);
				vertices[i].x = -1000;
				vertices[i].y = -1000;
				vertices[i].z = -1000;
			}
		}


		this.geometry.elementsNeedUpdate = true;
		// this.geometry.buffersNeedUpdate = true;


		this.geometry.computeFaceNormals();
		this.geometry.computeVertexNormals();

		this.buildNeighborMap();



	};

}
