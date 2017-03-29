Assignment 2 complex data

Patikas Nikolaos, 1747

22-4-16

Usage
Just compile in the folder with make. Compiler must support the c++11 flag.
run with: 
	$ ./shortest <source-id> <destination-id>

Implementation

Functions work as intended optimized for performance with the use of STL indexes. Additionally some there are provided some timing functions to give an cost per iteration idea. If you want to view the paths you 'll have to tweak the view_path variable.
	  functions:
		distance( Point p1 ,Point p2) :
			  gives the euclidean distance between the geographical locations between two nodes.

		distance(Location* source,
				   Location* dest,
				   unordered_map<Location*, unordered_map<Location*, double> >& landmarks):
			estimates distance with the use of landmarks.

		dijkstra(Location *source, Location *target):
			finds best path using the dijkstra method

		dijkstra(Location *source,
				  multimap<double, Location*>& nearest,
				  unordered_map<Location*, double>& dist)
			used for landmark node to find shortest path cost to all other nodes of the grid.

		aStarSearch(Location *start, Location *end)

		       finds the best path from start to end using A star search and the euclidean distance.

		findLandmarks(int N)
			gets N most distant landmarks and their arrays of distance. Works with random seed node v0

		oneAstarSearch(Location *start,
					Location *end,
					unordered_map<Location*, unordered_map<Location*, double> >& dist)
			a star search using heuristic based on landmarks.


Results

We see a domination of 1A* in time of execution for distant nodes and in iterations for all nodes. Although it has a fixed cost for building the landmark arrays this cost pays off for distant nodes. Additionally the landmark arrays are reusable. Second in rank comes A* star search with the euclidean distance which has smaller iteration cost than 1A* but it takes more iterations to find the shortest path, and last comes the dijkstra algorithm which is seems to be more suited for forwarding arrays as it has no apparent direction to the destination and it takes much more iterations to reach its destination than the other two.

	Thoughts:
		-seed node v0 for landmarks does not matter. Same landmarks will end up in list in most cases.
		-1A* search is the best search algorithm for nodes that are not very close.
		-Build 10 landmark arrays has 0.25 seconds overhead. Time is not included in 1A* algorithm execution time.

Examples:

Without paths:

Near nodes (15 ->205) winner A star

$ ./shortest 15 205
Reading nodes from file 'cal.cnode' and edges from file 'cal.cedge'
Importing  finished. 21048 locations in the grid
Initial seed node 11720
Landmark 1	 Node: 11720	Cost: 0
Landmark 2	 Node: 80	Cost: 17.4865
Landmark 3	 Node: 20600	Cost: 17.4865
Landmark 4	 Node: 482	Cost: 37.5272
Landmark 5	 Node: 20618	Cost: 47.4041
Landmark 6	 Node: 8264	Cost: 57.4654
Landmark 7	 Node: 16234	Cost: 62.5239
Landmark 8	 Node: 79	Cost: 77.1457
Landmark 9	 Node: 20599	Cost: 82.3949
Landmark 10	 Node: 483	Cost: 97.4692
Seconds to build landmark arrays: 0.253487
1A star Search ::	(15 -> 205)	Iterations 534
Cost: 3.3897 Seconds: 0.014071
A star Search ::	(15 -> 205)	Iterations 835
Cost: 3.3897 Seconds: 0.019981
Dijkstra Search ::	(15 -> 205)	Iterations 3221
Cost: 3.3897 Seconds: 0.06908

Far Nodes: Big Winner 1A*
$ ./shortest 15 20221
Reading nodes from file 'cal.cnode' and edges from file 'cal.cedge'
Importing  finished. 21048 locations in the grid
Initial seed node 1499
Landmark 1	 Node: 1499	Cost: 0
Landmark 2	 Node: 20600	Cost: 28.0125
Landmark 3	 Node: 80	Cost: 28.0125
Landmark 4	 Node: 16234	Cost: 34.1469
Landmark 5	 Node: 482	Cost: 50.2875
Landmark 6	 Node: 8264	Cost: 50.9725
Landmark 7	 Node: 20618	Cost: 63.5501
Landmark 8	 Node: 79	Cost: 77.9845
Landmark 9	 Node: 20599	Cost: 82.3949
Landmark 10	 Node: 483	Cost: 97.4692
Seconds to build landmark arrays: 0.255848
1A star Search ::	(15 -> 20221)	Iterations 950
Cost: 13.4859 Seconds: 0.076979
A star Search ::	(15 -> 20221)	Iterations 12002
Cost: 13.4859 Seconds: 0.844181
Dijkstra Search ::	(15 -> 20221)	Iterations 20630
Cost: 13.4859 Seconds: 1.12986



With paths:

Near nodes: (180 -> 20)
$ ./shortest 180 20
Reading nodes from file 'cal.cnode' and edges from file 'cal.cedge'
Importing  finished. 21048 locations in the grid
Initial seed node 5486
Landmark 1	 Node: 5486	Cost: 0
Landmark 2	 Node: 20600	Cost: 22.6423
Landmark 3	 Node: 80	Cost: 22.6423
Landmark 4	 Node: 16234	Cost: 34.1469
Landmark 5	 Node: 482	Cost: 50.2875
Landmark 6	 Node: 8264	Cost: 50.9725
Landmark 7	 Node: 20618	Cost: 63.5501
Landmark 8	 Node: 79	Cost: 77.9845
Landmark 9	 Node: 20599	Cost: 82.3949
Landmark 10	 Node: 483	Cost: 97.4692
Seconds to build landmark arrays: 0.253577
1A star Search ::	(180 -> 20)	Iterations 117
(180,	181,	182,	183,	184,	185,	186,	187,	188,	189,	190,	191,	192,	193,	194,	195,	196,	197,	198,	199,	200,	201,	202,	203,	83,	84,	85,	86,	87,	88,	89,	90,	91,	92,	93,	94,	95,	96,	97,	98,	82,	81,	39,	38,	37,	9,	10,	11,	12,	13,	14,	15,	16,	17,	18,	19,	20)
Cost: 0.865029 Seconds: 0.001203
A star Search ::	(180 -> 20)	Iterations 105
(180,	181,	182,	183,	184,	185,	186,	187,	188,	189,	190,	191,	192,	193,	194,	195,	196,	197,	198,	199,	200,	201,	202,	203,	83,	84,	85,	86,	87,	88,	89,	90,	91,	92,	93,	94,	95,	96,	97,	98,	82,	81,	39,	38,	37,	9,	10,	11,	12,	13,	14,	15,	16,	17,	18,	19,	20)
Cost: 0.865029 Seconds: 0.000692
Dijkstra Search ::	(180 -> 20)	Iterations 245
(180,	181,	182,	183,	184,	185,	186,	187,	188,	189,	190,	191,	192,	193,	194,	195,	196,	197,	198,	199,	200,	201,	202,	203,	83,	84,	85,	86,	87,	88,	89,	90,	91,	92,	93,	94,	95,	96,	97,	98,	82,	81,	39,	38,	37,	9,	10,	11,	12,	13,	14,	15,	16,	17,	18,	19,	20)
Cost: 0.865029 Seconds: 0.001429
