//============================================================================
// Name        : shortest.cpp
// Author      : Nikolaos Patikas 1747
//============================================================================

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <queue>
#include <functional>
#include <cmath>
#include <sys/time.h>
#include <fstream>
#include <set>
#include <unordered_set>
#include <utility>
#include <cfloat>
#define D 2
using namespace std;


// Variables for paths
bool dump_path = false;
bool view_path = true;

//-------------------------------
//-------Class Declarations------
//-------------------------------

class Point {
public:
    double x, y;
    Point(double x, double y) :
	x(x), y(y) {
    }
};

class Location {
public:
    int id;
    Point coordinates;
    unordered_map<Location*, double> neighbors;
    priority_queue<pair<Location*, double>, vector<pair<Location*, double> > > distance_matrix;

    Location(Point l, int id) :
	coordinates(l), id(id) {
    }
    void addEdge(Location* v, double d);
	
    //hash function for a location node
    size_t operator ()(const Location* l) {
	return hash<double>()(l->coordinates.x)
	    ^ hash<double>()(l->coordinates.y);
    }

};

void Location::addEdge(Location*v, double d) {
    neighbors.insert(make_pair(v, d));
}

// hash function for class Location
struct locationCompare {
    bool operator ()(const Location* l1, const Location* l2) const {
	return (l1->coordinates.x == l2->coordinates.x)
	    && (l1->coordinates.y == l2->coordinates.y);
    }
};

//-------------------------------
//---End of Class Declarations---
//-------------------------------

map<int, Location*> locations; // as exercise asked
vector<int> path;


//importing nodes
void readFiles(string nodes_file, string edges_file) {
    ifstream node_fp(nodes_file.c_str());

    cout << "Reading nodes from file \'" << nodes_file
	 << "\' and edges from file \'" << edges_file << "\'" << endl;

    if(!node_fp.is_open()){
	cout<<"File "<<nodes_file<<" not found"<<endl;
	exit(-1);
    }
    while (!node_fp.eof()) {
	int id;
	double longtitude, latitude;
	node_fp >> id >> longtitude >> latitude;
	Point p(longtitude, latitude);
	locations[id] = new Location(p, id);

    }

    ifstream edge_fp(edges_file.c_str());
    if(!edge_fp.is_open()){
	cout<<"File "<<edges_file<<" not found"<<endl;
	exit(-1);	
    }
    while (!edge_fp.eof()) {
	int edge_id, node1, node2;
	double distance;
	edge_fp >> edge_id >> node1 >> node2 >> distance;
	locations[node1]->addEdge(locations[node2], distance);
	locations[node2]->addEdge(locations[node1], distance);
    }
    cout << "Importing  finished. " << locations.size()
	 << " locations in the grid" << endl;

}

// Custom landmark distance function
double distance(Location* source, Location* dest, unordered_map<Location*, unordered_map<Location*, double> >& landmarks) {

    double max = 0;
    for(auto& it : landmarks){
	auto landmark = it.first;
	auto& dist = it.second;
	double h = fabs(dist[source] - dist[dest]);
	if(h > max)
	    max = h;
    }
    return max;

}

// Always underestimates the distance unless teleporting
double distance(Point p1, Point p2) {
    
    return sqrt(((p1.x - p2.x) * (p1.x - p2.x)) + ((p1.y - p2.y) * (p1.y - p2.y)));

}


void construct_path(Location*start, Location*node, unordered_map<Location*, pair<Location*, double> > came_from) {

    if (node != start) {
	construct_path(start, came_from[node].first, came_from);

    }
    if(dump_path){
	cout << "-> " << node->id;
    }
    path.push_back(node->id);

}


void construct_path(Location*start, Location*node, unordered_map<Location*, Location*> came_from) {

    if (node != start && node) {
	construct_path(start, came_from[node], came_from);

    }
    if(dump_path){
	cout << " -> " << node->id;
    }
    path.push_back(node->id);
}

//Traditional dijkstra from node to node 2nd exercise
double dijkstra(Location *source, Location *target) {
    multimap<double, Location*> Q;
    unordered_map<Location*, double> inv_Q;
    unordered_map<Location*, double> dist;
    unordered_map<Location*, Location*> prev;
    dist[source] = 0;
    prev[source] = NULL;

    Q.insert(make_pair(0, source));
    inv_Q[source] = 0;
    int iterations = 0;

    while (!Q.empty()) {
	iterations++;
	Location *u = Q.begin()->second;
	Q.erase(Q.begin());
	if (u == target) {
	    cout << "Dijkstra Search ::\t(" << source->id << " -> " << target->id <<  ")\tIterations " << iterations << endl;
	    construct_path(source,u,prev);
	    return dist[u];

	}

	for (auto it = u->neighbors.begin(); it != u->neighbors.end(); ++it) {
	    auto v = it->first;
	    auto cost = it->second;
	    double alt_route = dist[u] + cost;
	    if (dist.find(v) == dist.end() || dist[v] > alt_route) {
		dist[v] = alt_route;
		prev[v] = u;

		//decrease priority of node v
		auto node = inv_Q.find(v);
		if (node == inv_Q.end()) {
		    Q.insert(make_pair(alt_route, v));
		    inv_Q[v] = alt_route;
		} else {
		    auto old_cost = inv_Q[v];
		    for (auto it2 = Q.find(old_cost); it2 != Q.end(); ++it2) {
			if (it2->second == v) {
			    Q.erase(it2);
			    break;
			}

		    }
		    Q.insert(make_pair(alt_route, v));
		    inv_Q[v] = alt_route;

		}

	    }
	}
    }

}


// Dijkstra for landmarks
void dijkstra(Location *source, multimap<double, Location*>& nearest, unordered_map<Location*, double>& dist) {
    
    multimap<double, Location*> Q;
    unordered_map<Location*, double> inv_Q;
    unordered_map<Location*, Location*> prev;
    
    dist[source] = 0;
    prev[source] = NULL;
    
    Q.insert(make_pair(0, source));
    inv_Q[source] = 0;

    while (!Q.empty()) {

	Location *u = Q.begin()->second;
	Q.erase(Q.begin());

	for (auto it : u->neighbors) {
	    auto v = it.first;
	    auto cost = it.second;
	    auto alt_route = dist[u] + cost;
	    if (dist.find(v) == dist.end() || dist[v] > alt_route) {
		dist[v] = alt_route;
		prev[v] = u;

		//decrease priority of node v
		auto node = inv_Q.find(v);
		if (node == inv_Q.end()) {
		    Q.insert(make_pair(alt_route, v));
		    inv_Q[v] = alt_route;
		} else {
		    auto old_cost = inv_Q[v];
		    for (auto it2 = Q.find(old_cost); it2 != Q.end(); ++it2)
			if (it2->second == v) {
			    Q.erase(it2);
			    break;
			}
		    Q.insert(make_pair(alt_route, v));
		    inv_Q[v] = alt_route;

		}

	    }
	}
    }

    // give an ordered structure of nearest nodes according to distance
    for (auto it : dist)
	nearest.insert(make_pair(it.second, it.first));

}


double aStarSearch(Location *start, Location *end) {

    multimap<double, Location*> open_set;
    unordered_map<Location*, double> open_set_inverse;
    set<Location*> closed_set;
    unordered_map<Location*, pair<Location*, double> > came_from;
    int iterations = 0;

    open_set.insert(make_pair(distance(start->coordinates, end->coordinates), start));
    open_set_inverse.insert(make_pair(start, distance(start->coordinates, end->coordinates)));
    came_from.insert(make_pair(start, make_pair(start, 0)));
    
    while (!open_set.empty()) {
	Location *current = open_set.begin()->second;
	double gscore_current = came_from[current].second;
	
	iterations++;
	closed_set.insert(current);
	open_set.erase(open_set.begin());
	open_set_inverse.erase(current);

	if (current == end) {
	    cout << "A star Search ::\t(" << start->id << " -> " << end->id <<  ")\tIterations " << iterations << endl;;
	    construct_path(start, end, came_from);
	    return gscore_current;
	}
	
	for (auto &it : current->neighbors) {
	    
	    Location* neighbor = it.first;
	    double edge_cost = it.second;

	    if (closed_set.find(neighbor) != closed_set.end())
		continue;

	    double g_score = gscore_current + edge_cost;
	    auto node_in_open_set = open_set_inverse.find(neighbor);
	    auto current_score = came_from.find(neighbor);

	    if (node_in_open_set == open_set_inverse.end() || (current_score->second.second > g_score ) ) {
		//cout<<"I love programming shit"<<endl;
		double f_score = g_score + distance(neighbor->coordinates, end->coordinates);

		came_from.erase(neighbor);
	        came_from.insert(make_pair(neighbor, make_pair(current, g_score)));

		if (node_in_open_set != open_set_inverse.end()) {
		    for (auto it2 = open_set.find(node_in_open_set->second); it2 != open_set.end(); ++it2) {
			if (it2->second == neighbor) {
			    open_set.erase(it2);
			    break;
			}
		    }

		    open_set_inverse.erase(neighbor);

		}
		open_set.insert(make_pair(f_score, neighbor));
		open_set_inverse.insert(make_pair(neighbor, f_score));
	    }

	}
    }
    cout << "Failure! Could not find a path!" << endl;
    return -1;

}
double seconds_elapsed(struct timeval *start, struct timeval *end){
    
    double elapsedTime = 0;
    
    elapsedTime = (end->tv_sec - start->tv_sec);
    elapsedTime += (end->tv_usec - start->tv_usec) /1000000.0;

    return elapsedTime;
}

unordered_map<Location*, unordered_map<Location*, double> > findLandmarks(int N) {

    //Select Random Landmark
    int v0=rand()%locations.size();
    
    if (locations.find(v0) != locations.end()) {
	cout << "Initial seed node " << v0 << endl;
    } else {
	cout << "Random selected landmark " << v0 << " doesn't exist" << endl;
	exit(1);
    }
    
    //Final landmarks product
    unordered_map<Location*, unordered_map<Location*, double> > landmarks;
    //Indexed by distance array
    unordered_map<Location*, multimap<double, Location*> > ordered_distance;

    //Locations with reversed indexes
    
    unordered_map<Location*, double> total_dist;
    //Initialize distance counters
        
    //Find Other Landmarks
    Location *current = locations[v0];

    double max = 0;
    Location* max_node = NULL;
    
    for (int i = 0; i < N; i++) {

	unordered_map<Location*, double> dist;
	multimap<double, Location*> nearest;

	dijkstra(current, nearest, dist);
	landmarks[current] = dist;
	cout <<"Landmark "<< i + 1 << "\t Node: " << current->id << "\tCost: " << max << endl;
	if (i == 1) {

	    auto it = nearest.rbegin();
	    current = it->second;
	    for(auto& it : dist)
		total_dist[it.first] = it.second;
	    continue;

	}
	
	for (auto& it : dist)
	    total_dist[it.first] += it.second;
	     



	for (auto & it : dist){
	    double temp = total_dist[it.first] + it.second;
	    if (max < temp && landmarks.find(it.first) == landmarks.end()){
		max = temp;
		max_node = it.first;
	    }
	}
	current = max_node;
	
	

    }
    return landmarks;
}

 

double oneAstarSearch(Location *start, Location *end, unordered_map<Location*, unordered_map<Location*, double> >& dist) {

    multimap<double, Location*> open_set;
    unordered_map<Location*, double> open_set_inverse;
    set<Location*> closed_set;
    unordered_map<Location*, pair<Location*, double> > came_from;
    int iterations = 0;

    open_set.insert(make_pair(distance(start, end, dist), start));
    open_set_inverse.insert(make_pair(start, distance(start, end, dist)));
    came_from.insert(make_pair(start, make_pair(start, 0)));
    
    while (!open_set.empty()) {
	Location *current = open_set.begin()->second;
	double gscore_current = came_from[current].second;
	
	iterations++;
	closed_set.insert(current);
	open_set.erase(open_set.begin());
	open_set_inverse.erase(current);

	if (current == end) {
	    cout << "1A star Search ::\t(" << start->id << " -> " << end->id <<  ")\tIterations " << iterations << endl;
	    construct_path(start, end, came_from);
	    return gscore_current;
	}
	
	for (auto &it : current->neighbors) {
	    
	    Location* neighbor = it.first;
	    double edge_cost = it.second;

	    if (closed_set.find(neighbor) != closed_set.end())
		continue;

	    double g_score = gscore_current + edge_cost;
	    auto node_in_open_set = open_set_inverse.find(neighbor);
	    auto current_score = came_from.find(neighbor);

	    if (node_in_open_set == open_set_inverse.end() || (current_score->second.second > g_score ) ) {
		//cout<<"I love programming shit"<<endl;
		double f_score = g_score + distance(neighbor, end, dist);

		came_from.erase(neighbor);
	        came_from.insert(make_pair(neighbor, make_pair(current, g_score)));

		if (node_in_open_set != open_set_inverse.end()) {
		    for (auto it2 = open_set.find(node_in_open_set->second); it2 != open_set.end(); ++it2) {
			if (it2->second == neighbor) {
			    open_set.erase(it2);
			    break;
			}
		    }

		    open_set_inverse.erase(neighbor);

		}
		open_set.insert(make_pair(f_score, neighbor));
		open_set_inverse.insert(make_pair(neighbor, f_score));
	    }

	}
    }
    cout << "Failure! Could not find a path!" << endl;
    return -1;

}

void dumpPath(vector <int>& p){
    if(!view_path)
	return;
    cout << "(";

    for(vector<int>::const_iterator i = p.begin(); i != p.end(); ++i){
	cout << *i;
	if(p.end() != i + 1 ){
	    cout <<",\t";
	    }
    }

    cout << ")" << endl;


}


void findRoute(unordered_map<Location*, unordered_map<Location*, double> >& landmarks , int source, int destination){

    double cost=0;
    struct timeval t1, t0;
    
    if (locations.find(source) == locations.end()) {
	cout << "Source node " << source << " doesn't exist" << endl;
	exit(1);
    }
    if (locations.find(destination) == locations.end()) {
	cout << "Destination node " << destination << " doesn't exist" << endl;
	exit(1);
    }
    gettimeofday(&t0, NULL);
    cost = oneAstarSearch(locations[source],locations[destination],landmarks);
    auto onestarpath = path;
    path.clear();
    gettimeofday(&t1, NULL);
    dumpPath(onestarpath);
    cout << "Cost: " << cost <<" Seconds: " << seconds_elapsed(&t0, &t1) <<endl;
    
    gettimeofday(&t0, NULL);
    cost = aStarSearch(locations[source], locations[destination]);
    auto astarpath = path;
    path.clear();
    gettimeofday(&t1, NULL);
    dumpPath(onestarpath);
    cout << "Cost: " << cost <<" Seconds: " << seconds_elapsed(&t0, &t1) <<endl;
    
    gettimeofday(&t0, NULL);
    cost = dijkstra(locations[source], locations[destination]);
    auto dijkstrapath = path;
    path.clear();
    gettimeofday(&t1, NULL);
    dumpPath(dijkstrapath);
    cout << "Cost: " << cost <<" Seconds: " << seconds_elapsed(&t0, &t1) <<endl;


}

int main(int argc,char*argv[]){

    map<double, double> distance_matrix;
    char nfile[]="cal.cnode";
    char efile[]="cal.cedge";
    vector<int> nodes;
    srand(time(NULL));
    if(argc!=3){
	cout<<"Default files:"<<nfile<<" "<<efile<<endl;
	cout<<"Usage: "<<argv[0]<<" <source-id> <dest-id>"<<endl;
    }else{
	readFiles(nfile,efile);
	struct timeval t1, t0;
	gettimeofday(&t0, NULL);
	auto landmarks = findLandmarks(10);
	gettimeofday(&t1, NULL);
	cout << "Seconds to build landmark arrays: " << seconds_elapsed(&t0, &t1) <<endl;

	findRoute(landmarks, atoi(argv[1]), atoi(argv[2]));
    }

    return 0;

}



