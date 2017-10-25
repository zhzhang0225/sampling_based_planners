#ifndef PRMGRAPH_H
#define PRMGRAPH_H

#include <math.h>
#include <float.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>

/**************************************************** 
Class for creating and maintaining a graph structure 
as the representation of Probabilistic Roadmap
****************************************************/

#define PI 3.141592654

using namespace std;

class PRMGraph {
private:
	//node_id (initialized to zero!)
	int node_id_;
	//record number of degree of freedom
	int num_dof_;
	//represent vertices as array of joint states 
	vector<vector<double>> vertices_;
	//represent edges as adjacency list 
	//1st dim = vertex, 2nd dim = neighbors of the vertex
	vector<vector<int>> edges_;

public:
	//default constructor
	PRMGraph(int numofDOFs);

	int getCurrentNodeID();

	void addVertex(vector<double>& vertex);

	vector<double> getNodeConfig(int node_id);

	vector<int> getNeighborsID(int node_id);

	void addEdge(int parent_id, int child_id);

	inline double calculateDistance(vector<double>& v1, vector<double>& v2)
	{
		//calculate L2 (norm) distance between two configurations (vertices)
		double total_dist = 0.0;
		double joint_dist = 0.0;
		for (int i=0; i<(this->num_dof_); i++) {
			joint_dist = fabs(v1[i]-v2[i]);
			total_dist += (joint_dist < PI) ? pow(joint_dist,2) : pow((2*PI - joint_dist),2);
		}
		return pow(total_dist, 0.5);
	};

	vector<int> findKNN(vector<double>& new_vertex, int k);

	int getNearestVertex(vector<double>& new_vertex);
};

#endif