#include "prm_graph.h"

PRMGraph::PRMGraph(int numofDOFs)
: num_dof_(numofDOFs),
  node_id_(0) { }

int PRMGraph::getCurrentNodeID() 
{
	return this->node_id_;
}

void PRMGraph::addVertex(vector<double>& vertex)
{
	// add one row
	(this->vertices_).push_back(vector<double>());
	// fill in the new row 
	vertices_[this->node_id_] = vertex;

	// add place-holder  
	(this->edges_).push_back(vector<int>());

	(this->node_id_)++;
}

vector<double> PRMGraph::getNodeConfig(int node_id)
{
	return (this->vertices_)[node_id];
}

vector<int> PRMGraph::getNeighborsID(int node_id)
{
	return (this->edges_)[node_id];
}

void PRMGraph::addEdge(int vertex1_id, int vertex2_id)
{
	(this->edges_)[vertex1_id].push_back(vertex2_id);
	(this->edges_)[vertex2_id].push_back(vertex1_id);
	//std::cout << (this->edges_)[vertex1_id].size() << " " << (this->edges_)[vertex2_id].size() << std::endl;
}

vector<int> PRMGraph::findKNN(vector<double>& new_vertex, int k)
{
	vector<pair<int, double>> dist;
	vector<int> knn_id;
	// perform linear search to get K nearest neighbors 

	// calculate distance for all vertices
	for (int i=0; i<(this->vertices_).size(); i++) {
		auto curr_dist = make_pair(i, this->calculateDistance(new_vertex, (this->vertices_)[i]));
		dist.push_back(curr_dist);
	}
	// sort the distance 
	sort(dist.begin(), dist.end(),
		[](const pair<int, double>& p1, const pair<int, double>& p2)
		{
			return (p1.second < p2.second);
		});

	// find id of k nearest neighbors
	for (int i=0; i<dist.size(); i++) {
		if (dist[i].second > 0.1) knn_id.push_back(dist[i].first);
		if (knn_id.size() >= k) break;
	}
	return knn_id;
}

int PRMGraph::getNearestVertex(vector<double>& new_vertex) 
{
	double min_dist = DBL_MAX;
	int min_index = -1;
	double curr_dist = 0.0;
	// perform linear search to get the nearest neighbor
	for (int i=0; i<(this->vertices_).size(); i++) {
		curr_dist = this->calculateDistance(new_vertex, (this->vertices_)[i]);
		if (curr_dist < min_dist) {
			min_dist = curr_dist;
			min_index = i;
		}
	}
	return min_index;
}