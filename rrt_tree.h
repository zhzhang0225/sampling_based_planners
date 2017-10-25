#ifndef RRTTREE_H
#define RRTTREE_H

#include <math.h>
#include <float.h>
#include <vector>
#include <unordered_map>
#include <assert.h>

/**************************************************** 
Class for creating and maintaining a growing tree 
structure used by RRT/RRT-connect/RRT* algorithms 
****************************************************/

#define PI 3.141592654

using namespace std;

class RRTTree {
private:
	//node_id (initialized to zero!)
	int node_id_;
	//record number of degree of freedom
	int num_dof_;
	//represent vertices as array of joint states 
	vector<vector<double>> vertices_;
	//represent edge by its vertices' id
	//key = child_id, value = parent_id
	unordered_map<int, int> edges_;
	//cost associated with each vertex (for RRT* only)
	vector<double> costs_;

public:
	//default constructor 
	RRTTree(int numofDOFs);

	int getNodeID();

	vector<double> getNodeConfig(int node_id);

	//int getParentVertex(int node_id);

	void removeEdge(int child_node_id);

	void addVertex(vector<double>& vertex);

	double getVertexCost(int node_id);

	void setVertexCost(int node_id, double cost);

	void addEdge(int parent_id, int child_id);

	inline double calculateDistance(vector<double>& v1, vector<double>& v2)
	{
		//calculate distance between two configurations (vertices)
		//as the sum of differences between all joint angles
		double total_dist = 0.0;
		double joint_dist = 0.0;
		for (int i=0; i<(this->num_dof_); i++) {
			joint_dist = fabs(v1[i]-v2[i]);
			total_dist += (joint_dist < PI) ? joint_dist : (2*PI - joint_dist);
		}
		return total_dist;
	};

	int getNearestVertex(vector<double>& new_vertex);

	vector<int> getNearVertices(int node_id, double radius);

	vector<int> returnPlan();

	vector<int> returnPlan(int node_id);

};

#endif 