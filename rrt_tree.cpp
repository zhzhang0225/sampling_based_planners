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

void RRTTree::addVertex(vector<double>& vertex)
{
	assert(vertex.size() == this->num_dof_);

	// add one row
	(this->vertices_).push_back(vector<double>());
	// fill in the new row 
	vertices_[this->node_id_] = vertex;

	(this->node_id_)++;
}

void RRTTree::addEdge(int parent_id, int child_id)
{
	assert((this->edges_).find(child_id) == (this->edges_).end());
	(this->edges_)[child_id] = parent_id;
}

// inline double RRTTree::calculateDistance(vector<double>& v1, vector<double>& v2)
// {
// 	//calculate distance between two configurations (vertices)
// 	//as the sum of differences between all joint angles
// 	double total_dist = 0.0;
// 	double joint_dist = 0.0;
// 	for (int i=0; i<(this->num_dof_); i++) {
// 		joint_dist = fabs(v1[i]-v2[i]);
// 		total_dist += (joint_dist < PI) ? joint_dist : (2*PI - joint_dist);
// 	}
// 	return total_dist;
// }

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