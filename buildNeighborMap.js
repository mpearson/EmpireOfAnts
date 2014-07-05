
var pairs = [
	'D', 'C',
	'A', 'B',
	'E', 'F',
	'F', 'A',
	'E', 'D',
	'C', 'B'
];



var buildCycle = function(pairs) {
	var count = pairs.length/2;
	var cycle = Array(count);

	cycle[0] = pairs[0];
	cycle[1] = pairs[1];

	var position = 1;

	// search for a pair containing the current vertex, and add its comrade to the end of the cycle
	while(position < count - 1) {
		// we can skip the first 2 since they've already been added
		for(var j=2; j<pairs.length; j+=2) {
			if(pairs[j] === cycle[position])
				cycle[++position] = pairs[j+1];
			else if(pairs[j+1] === cycle[position])
				cycle[++position] = pairs[j];
			else
				continue;
			pairs[j] = pairs[j+1] = null;	// once it's added to the cycle, we don't want to find it again
			break;
		}
	}
	return cycle;
};

cycle = buildCycle(pairs);
