#include "rrt_tree.h"

RRTTree::RRTTree(int numofDOFs) 
: num_dof_(numofDOFs),
  node_id_(0) { }


int RRTTree::getNodeID() 
{
	return this->node_id_;
}

vector<double> RRTTree::getNodeConfig(int node_id)
{
	return (this->vertices_)[node_id];
}

// int RRTTree::getParentVertex(int node_id)
// {
// 	return (this->edges_)[node_id];
// }

void RRTTree::removeEdge(int child_node_id)
{
	(this->edges_).erase(child_node_id);
	return;
}

void RRTTree::addVertex(vector<double>& vertex)
{
	assert(vertex.size() == this->num_dof_);

	// add one row
	(this->vertices_).push_back(vector<double>());
	// fill in the new row 
	vertices_[this->node_id_] = vertex;

	(this->node_id_)++;
}

double RRTTree::getVertexCost(int node_id)
{
	assert(node_id < (this.costs_).size());
	return (this->costs_)[node_id];
}

void RRTTree::setVertexCost(int node_id, double cost) 
{
	if (node_id < (this->costs_).size()) {
		// update cost for exisiting vertex
		(this->costs_)[node_id] = cost;
	} else {
		// record cost for new vertex
		assert(node_id == this->node_id_-1);
		(this->costs_).push_back(cost);
	}
}

void RRTTree::addEdge(int parent_id, int child_id)
{
	assert((this->edges_).find(child_id) == (this->edges_).end());
	(this->edges_)[child_id] = parent_id;
	//cout << (this->edges_).size() << ":  " << parent_id << ", " << child_id << endl;
}

int RRTTree::getNearestVertex(vector<double>& new_vertex) 
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

vector<int> RRTTree::getNearVertices(int node_id, double radius)
{
	//perform linear search to get neighbors within distance r to current vertex
	vector<int> neighbors_id;
	double curr_dist = 0.0;
	for (int i=0; i<(this->vertices_).size(); i++) {
		curr_dist = this->calculateDistance((this->vertices_)[node_id], (this->vertices_)[i]);
		if (curr_dist <= radius) neighbors_id.push_back(i);
	}
	return neighbors_id;
}

vector<int> RRTTree::returnPlan()
{
	// return the vertices' ids along the path
	vector<int> plan_vertices;
	int index_iter = this->node_id_-1;
	while (index_iter) {
		plan_vertices.push_back(index_iter);
		index_iter = (this->edges_)[index_iter];
	}
	// add the start configuration
	plan_vertices.push_back(0);
	return plan_vertices;
}

vector<int> RRTTree::returnPlan(int node_id)
{
	// return the vertices' ids along the path
	vector<int> plan_vertices;
	int index_iter = node_id;
	while (index_iter) {
		plan_vertices.push_back(index_iter);
		index_iter = (this->edges_)[index_iter];
	}
	// add the start configuration
	plan_vertices.push_back(0);
	return plan_vertices;
}