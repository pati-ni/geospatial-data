//============================================================================
// Name        : MYE014Assignment2.cpp
// Author      : Nikos Patikas 1747
//============================================================================

#include <iostream>
#include <vector>
#include <sstream>
#include <map>
#include <unordered_map>
#include <queue>
#include <functional>
#include <cmath>
#include <fstream>
#include <set>
#include <unordered_set>
#include <utility>
#include <cfloat>
#include <sys/time.h>
#define D 2
#define INF DBL_MAX
#define NRA_ALG 2
#define TA_ALG 1
#define SUM 1
#define MAXIMUM 2
using namespace std;

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
		double dist;
		edge_fp >> edge_id >> node1 >> node2 >> dist;
		locations[node1]->addEdge(locations[node2], dist);
		locations[node2]->addEdge(locations[node1], dist);
	}
	cout << "Importing  finished. " << locations.size()
			<< " locations in the grid" << endl;

}
// Always underestimates the distance
double distance(Point p1, Point p2) {
	return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

bool sum = true;

void construct_path(Location*start, Location*node,
		unordered_map<Location*, pair<Location*, double> > came_from) {

	if (node != start) {
		construct_path(start, came_from[node].first, came_from);

	}
	cout << "-> " << node->id;

}

void dijkstra(Location *source, multimap<double, Location*>& nearest) {
	multimap<double, Location*> Q;
	unordered_map<Location*, double> inv_Q;
	unordered_map<Location*, double> dist;
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

	// give an ordered structure of the
	for (auto it : dist) {
		nearest.insert(make_pair(it.second, it.first));

	}

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
	    //cout << "A star Search ::\t(" << start->id << " -> " << end->id <<  ")\tIterations " << iterations << endl;;
	    //construct_path(start, end, came_from);
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

void dijkstra(Location *source, multimap<double, Location*>& nearest,
		unordered_map<Location*, double>& dist) {
	multimap<double, Location*> Q;
	unordered_map<Location*, double> inv_Q;
	//unordered_map<Location*,double> dist;
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

	// give an ordered structure of the
	for (auto it : dist) {
		nearest.insert(make_pair(it.second, it.first));

	}

}


void TA(vector<int>& nodes_id, bool sum) {

	double best_gamma = INF;
	Location* best_node = NULL;
	unordered_map<Location*,multimap<double, Location*> > nearest;
	unordered_map <Location* , multimap< double , Location* >::iterator > iterators;
	unordered_map <Location* , unordered_map< Location* , double > > dist;
	for (auto id : nodes_id) {
		multimap<double, Location*> near;
		unordered_map< Location* , double > c_dist;
		dijkstra(locations[id], near, c_dist);
		nearest.insert(make_pair(locations[id],near));
		iterators.insert(make_pair(locations[id],nearest[locations[id]].begin()));
		dist[locations[id]] = c_dist;
	}
	int iterations=0;
	while(1){
		//cout<<iterations<<endl;
		bool end_flag = true;
		//round robin
		for (auto& it : nearest) {
			auto& current_iterator = iterators[it.first];
			Location* near = current_iterator->second;
			double gamma = current_iterator->first;

			if(current_iterator==nearest[it.first].end() || gamma >= best_gamma){
				continue;
			}
			end_flag = false;
			iterations++;
			current_iterator++;

			for (auto id2 : nodes_id) {
				if (it.first->id == id2){
					continue;
				}
				//double cost = aStarSearch(locations[id2], near);

				double cost = dist[locations[id2]][near];
				if (sum){
					gamma += cost;
				}else{
					gamma = cost > gamma ? cost : gamma;
				}

			}
			if (gamma < best_gamma) {

				best_gamma = gamma;
				best_node = near;
				//cout<<gamma<<" "<<near->id<<endl;
			}
		}
		if(end_flag){
			break;
		}
	}

	double avg_cost=0;
	if(!sum){
		for(auto& n:nodes_id){
			avg_cost+=aStarSearch(locations[n], best_node);
		}
		avg_cost= avg_cost/nodes_id.size();
	}else{
		avg_cost = best_gamma/nodes_id.size();
	}
	for(auto id :nodes_id){
		cout << id <<" -> "<< best_node->id <<" Cost "<< dist[locations[id]][best_node]<<endl;

	}
	cout <<"Best meeting point is node :#"<< best_node->id <<" with average cost "<<avg_cost <<endl;
	cout<< iterations<<" Iterations"<<endl;
}


void NRA(vector<int> nodes_id, bool sum) {

	vector<multimap<double, Location*> > nearest_ordered;
	vector<unordered_map<Location*, double> > dist;

	unordered_map<Location*, unsigned int> counters;
	unordered_map<Location*, double> gamma;
	unordered_map<Location*, bool> visited;

	for (auto id : nodes_id) {
		multimap<double, Location*> t_or;
		unordered_map<Location*, double> t_d;
		dijkstra(locations[id], t_or, t_d);
		nearest_ordered.push_back(t_or);
		dist.push_back(t_d);
	}

	double best_gamma = INF;
	Location* vantage_point = NULL;
	int iterations = 0;
	while (1) {

		for (auto& user : nearest_ordered) {
			auto t = user.begin();
			double cost = t->first;
			Location* node = t->second;
			user.erase(t);
			//cout<<node->id<<endl;
			if (visited.find(node) != visited.end()) {
				continue;
			}
			iterations++;
			if (gamma.find(node) == gamma.end()) {
				gamma[node] = cost;
				counters[node] = 1;
			} else {
				counters[node] += 1;
				if (sum) {
					gamma[node] += cost;
				} else {
					gamma[node] = cost > gamma[node] ? cost : gamma[node];
				}
			}
			if (best_gamma > gamma[node] && counters[node] == nodes_id.size()) {
				best_gamma = gamma[node];
				vantage_point = node;
			}
			if (best_gamma <= gamma[node] || counters[node] == nodes_id.size()) {
				visited[node] = true;
				auto d = gamma.find(node);
				gamma.erase(d);
				auto c = counters.find(node);
				counters.erase(c);
			}

		}
		if (gamma.empty()) {
			break;
		}
		//cout<<vantage_point->id<<endl;
	}
	double avg_cost=0;
	if(!sum){
		for(auto& n:dist){
			avg_cost+=n[vantage_point];
		}
		avg_cost= avg_cost/nodes_id.size();
	}else{
		avg_cost = best_gamma/nodes_id.size();
	}
	for(int i = 0; i < nodes_id.size(); i++){
		cout << nodes_id[i] <<" -> "<< vantage_point->id <<" Cost "<< dist[i][vantage_point]<<endl;
	}
	cout <<"Best meeting point is node :#"<< vantage_point->id <<" with average cost "<<avg_cost <<endl;
	cout<< iterations<<" Iterations"<<endl;
}

void dijkstra(Location *source, Location *target) {
	multimap<double, Location*> Q;
	unordered_map<Location*, double> inv_Q;
	unordered_map<Location*, double> dist;
	unordered_map<Location*, Location*> prev;
	dist[source] = 0;
	prev[source] = NULL;

	Q.insert(make_pair(0, source));
	inv_Q[source] = 0;

	while (!Q.empty()) {

		Location *u = Q.begin()->second;
		Q.erase(Q.begin());
		if (u == target) {
			auto path = prev[u];
			while (path) {
				cout << "->" << path->id << " ";
				path = prev[path];
			}

		}

		for (auto it = u->neighbors.begin(); it != u->neighbors.end(); ++it) {
			auto v = it->first;
			auto cost = it->second;
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

void findSortestPath(int id, int id2) {
	if (locations.find(id) == locations.end()) {
		cout << "Node " << id << " not found" << endl;
		return;
	}
	if (locations.find(id2) == locations.end()) {
		cout << "Node " << id2 << " not found" << endl;
		return;
	}
	cout << "Node " << id << " to node " << id2 << ": " << endl;
	aStarSearch(locations[id], locations[id2]);

}

int invalid_arguments(string exec){
	cout<<"Usage: "<<exec<<" <node-1> <node-2> ... <node-n> <algorithm> <gamma-function>"<<endl;
	cout<<"Algorithm: "<<TA_ALG<<" = TA"<<NRA_ALG<<" = NRA "<<endl;
	cout<<"Gamma-function "<<SUM<<" = SUM"<<MAXIMUM<<" = MAX "<<endl;
	return -1;
}

double seconds_elapsed(struct timeval *start, struct timeval *end){

    double elapsedTime = 0;

    elapsedTime = (end->tv_sec - start->tv_sec);
    elapsedTime += (end->tv_usec - start->tv_usec) /1000000.0;

    return elapsedTime;
}
int main(int argc,char* argv[]) {

	map<double, double> distance_matrix;
	int algorithm_type = 0, gamma_function;
	vector<int>nodes;
	readFiles("cal.cnode", "cal.cedge");

	if(argc<4){
		return invalid_arguments("Νot enough args!");
	}else{
		for ( int i = 1 ; i < argc-2;i++){
			string s(argv[i]);
			stringstream ss(s);
			int id;
			ss>>id;
			if(ss.fail()){
				cout<<argv[i]<<" is not a valid node. Enter number!"<<endl;
				exit(-1);
			}
			nodes.push_back(id);
		}
		string s(argv[argc-2]);
		stringstream ss(s);
		ss>>algorithm_type;
		if(ss.fail()){
			return invalid_arguments(argv[argc-2]);

		}
		string s2(argv[argc-1]);
		stringstream ss2(s2);
		ss2>>gamma_function;
		if(ss.fail()){
			return invalid_arguments(argv[argc-1]);
		}
	}
	cout<<"Algorithm:";
	if(algorithm_type==NRA_ALG){
		cout<<"NRA";
	}else if (algorithm_type==TA_ALG){
		cout<<"TA";
	}
	cout<<" , Gamma function: ";
	if(gamma_function==SUM){
		cout<<"SUM";
	}else if(gamma_function==MAXIMUM){
		cout<<"MAXIMUM ";

	}
	cout<<endl<<"Node set:{ ";
	for(auto& n: nodes){
		cout<<n<<" ";
	}
	cout<<"}"<<endl;
	struct timeval t1, t0;
	gettimeofday(&t0, NULL);
	if(algorithm_type==NRA_ALG){
		NRA(nodes,gamma_function==SUM);
	}else if (algorithm_type==TA_ALG){
		TA(nodes,gamma_function==SUM);
	}
	gettimeofday(&t1, NULL);
	cout << "Seconds Elapsed: " << seconds_elapsed(&t0, &t1) <<endl;
	return 0;

}
