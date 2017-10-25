/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

#include <random>
#include <iostream>
#include <math.h>
#include <vector>
#include <assert.h>
#include <time.h>
#include <stack>
#include <unordered_map>
#include "rrt_tree.h"
#include "prm_graph.h"

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTCONNECT  1
#define RRTSTAR     2
#define PRM         3

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10

typedef struct {
  int X1, Y1;
  int X2, Y2;
  int Increment;
  int UsingYIndex;
  int DeltaX, DeltaY;
  int DTerm;
  int IncrE, IncrNE;
  int XIndex, YIndex;
  int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size)
{
    double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params)
{
  params->UsingYIndex = 0;

  if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
    (params->UsingYIndex)++;

  if (params->UsingYIndex)
    {
      params->Y1=p1x;
      params->X1=p1y;
      params->Y2=p2x;
      params->X2=p2y;
    }
  else
    {
      params->X1=p1x;
      params->Y1=p1y;
      params->X2=p2x;
      params->Y2=p2y;
    }

   if ((p2x - p1x) * (p2y - p1y) < 0)
    {
      params->Flipped = 1;
      params->Y1 = -params->Y1;
      params->Y2 = -params->Y2;
    }
  else
    params->Flipped = 0;

  if (params->X2 > params->X1)
    params->Increment = 1;
  else
    params->Increment = -1;

  params->DeltaX=params->X2-params->X1;
  params->DeltaY=params->Y2-params->Y1;

  params->IncrE=2*params->DeltaY*params->Increment;
  params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
  params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

  params->XIndex = params->X1;
  params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y)
{
  if (params->UsingYIndex)
    {
      *y = params->XIndex;
      *x = params->YIndex;
      if (params->Flipped)
        *x = -*x;
    }
  else
    {
      *x = params->XIndex;
      *y = params->YIndex;
      if (params->Flipped)
        *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params)
{
  if (params->XIndex == params->X2)
    {
      return 0;
    }
  params->XIndex += params->Increment;
  if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
    params->DTerm += params->IncrE;
  else
    {
      params->DTerm += params->IncrNE;
      params->YIndex += params->Increment;
    }
  return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
		   int x_size,
 		   int y_size)

{
	bresenham_param_t params;
	int nX, nY; 
    short unsigned int nX0, nY0, nX1, nY1;

    //printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
    
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

    //printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
            return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
		   int x_size, int y_size)
{
    double x0,y0,x1,y1;
    int i;
    
 	//iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
    y1 = 0;
	for(i = 0; i < numofDOFs; i++)
	{
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
				return 0;
	} 
	return 1;   
}

void incrementToJointAngles(int num_of_dof, double step_size, double** curr_ptr, std::vector<double>& goal)
{
	for (int i = 0; i < num_of_dof; i++) {
		double temp_dist = goal[i] - (*curr_ptr)[i];
		if (temp_dist > 0.0) {
			if (temp_dist < PI) (*curr_ptr)[i] += step_size;
			// remember to wrap the angle between 0 and 2pi
			else (*curr_ptr)[i] = (((*curr_ptr)[i] - step_size) > 0.0) ? ((*curr_ptr)[i] - step_size) 
				                                                : (2*PI + (*curr_ptr)[i] - step_size);
		}
		else {
			if (temp_dist > -PI) (*curr_ptr)[i] -= step_size; 
			// remember to wrap the angle between 0 and 2pi
			else (*curr_ptr)[i] = (((*curr_ptr)[i] + step_size) < 2*PI) ? ((*curr_ptr)[i] + step_size) 
				                                                 : ((*curr_ptr)[i] + step_size - 2*PI);
		}
	//assert((*curr_ptr)[i] > 0.0 && (*curr_ptr)[i] < 2*PI);
	}
	return;
}

static void planner(
		   double*	map,
		   int x_size,
 		   int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
	   int numofDOFs,
	   double*** plan,
	   int* planlength)
{
	//no plan by default
	*plan = NULL;
	*planlength = 0;
    
    //for now just do straight interpolation between start and goal checking for the validity of samples

    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("the arm is already at the goal\n");
        return;
    }
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    int firstinvalidconf = 1;
    for (i = 0; i < numofsamples; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf)
        {
            firstinvalidconf = 1;
            printf("ERROR: Invalid arm configuration!!!\n");
        }
    }    
    *planlength = numofsamples;
    
    return;
}


double calculateCost(std::vector<double>& config1, std::vector<double>& config2)
{
	// calcualte L1 cost between two configurations
	assert(config1.size() == config2.size());
	double total_cost = 0.0;
	for (int i = 0; i < config1.size(); i++) {
		total_cost += (fabs(config1[i]-config2[i])>PI) ? (2*PI-fabs(config1[i]-config2[i])) : fabs(config1[i]-config2[i]);
	}
	return total_cost;
}

static void plannerRRT(double*	map, 
					int x_size,
					int y_size,
					double* armstart_anglesV_rad,
					double* armgoal_anglesV_rad,
					int numofDOFs,
					double*** plan,
					int* planlength)
{
	std::uniform_real_distribution<double> goal_bias_generation(0.0, 1.0);
	std::uniform_real_distribution<double> joint_ang_generation(0.0, 2*PI);
	std::random_device rd;
	std::default_random_engine generator(rd());

	int epsilon = 60; //epsilon here specifies max. num of steps for each iteration
	double step_size = PI/90;

	int num_vertices = 0;

	// Create tree object 
	RRTTree tree(numofDOFs);

	// Add initial pose as vertex 
	std::vector<double> start_config(armstart_anglesV_rad, armstart_anglesV_rad+numofDOFs);
	tree.addVertex(start_config);
	tree.setVertexCost(tree.getNodeID()-1, 0.0);
	num_vertices++;

	// Create while loop until goal configuration is reached 
	bool goal_reached = false;
	
	while (!goal_reached) {
		// Initialize new/ sampled configuration
		std::vector<double> new_config(numofDOFs);

		// Introduce goal bias here
		double bias = goal_bias_generation(generator);
		if (bias > 0.80) {	
			// Set new_config to goal configuration
			copy(armgoal_anglesV_rad, armgoal_anglesV_rad+numofDOFs, new_config.begin());
		} 
		else {
			// Sample new arm configuration 
			//new_config[0] = joint_ang_generation(generator)/2;
			for (int i = 0; i < numofDOFs; i++)
				new_config[i] = joint_ang_generation(generator);
		}

		// Find nearest neighbor 
		int nn_index = tree.getNearestVertex(new_config);
		std::vector<double> nn_config = tree.getNodeConfig(nn_index);

		if (bias < 0.80 && tree.calculateDistance(new_config, nn_config) < 0.10) continue;

		// Attempt to extend to new sample 
		bool extend = true;
		double* temp_config = &nn_config[0];
		double* backup_config = (double*) malloc(numofDOFs*sizeof(double));
		int step_count = 0;

		while (extend) {
			// save configuration before increment 
			copy(temp_config, temp_config+numofDOFs, backup_config);
			// increment on each joint angle
			incrementToJointAngles(numofDOFs, step_size, &temp_config, new_config);

			// check collision 
			if (!IsValidArmConfiguration(temp_config, numofDOFs, map, x_size, y_size)) {
				copy(backup_config, backup_config+numofDOFs, temp_config);
				extend = false;
				break;
			} 
			step_count++;
			if (step_count > epsilon) extend = false;
 
			//vector<double> temp_config_vec(temp_config, temp_config+numofDOFs);
			//else if (tree.calculateDistance(temp_config_vec, new_config) < 0.1) extend = false;
		}

		free(backup_config);

		if (step_count > 0) {
			// Update ACUTUAL new configuration
			copy(temp_config, temp_config+numofDOFs, new_config.begin());
			// Add new_config to the tree
			tree.addVertex(new_config);
			tree.addEdge(nn_index, tree.getNodeID()-1);
			// Set cost for the new vertex
			nn_config = tree.getNodeConfig(nn_index);
			double cost_new = tree.getVertexCost(nn_index) + calculateCost(nn_config, new_config);
			tree.setVertexCost(tree.getNodeID()-1, cost_new);

			num_vertices++;
			if (num_vertices % 1000 == 0) std::cout << "Num of Vertices: " << num_vertices << std::endl;

			double dist_to_goal = 0.0;
			for (int i = 0; i < numofDOFs; i++) {
				double joint_dist = fabs(new_config[i] - armgoal_anglesV_rad[i]);
				dist_to_goal += (joint_dist > PI) ? (2*PI-joint_dist) : joint_dist;
			}
			if (dist_to_goal < 0.1) goal_reached = true; 
		}

		if (num_vertices > 50000) {
			std::cout << "Aborted due to timeout ... Please re-plan" << std::endl;
			break; 
		}
	}

    if (num_vertices <= 50000) {
	    std::cout << "generating the plan" << std::endl;
		//back-track the path and return the plan
		std::vector<int> plan_ids = tree.returnPlan();
		int path_length = plan_ids.size();

		*plan = (double**) malloc(path_length*sizeof(double*));
		for (int i = 0; i < path_length; i++) {
			(*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
			std::vector<double> config = tree.getNodeConfig(plan_ids[path_length-i-1]);
			//(*plan)[i] = &config[0]; is WRONG

			// for(int j = 0; j < numofDOFs; j++) std::cout << config[j] << "  ";
			// std::cout << std::endl;
			copy(config.begin(), config.end(), (*plan)[i]);
		}
		std::cout << "Cost to goal: " << tree.getVertexCost(tree.getNodeID()-1) << std::endl;
		*planlength = path_length;
	}

	return;
}


static void plannerRRTConnect(double* map, 
							  int x_size,
							  int y_size,
							  double* armstart_anglesV_rad,
							  double* armgoal_anglesV_rad,
							  int numofDOFs,
							  double*** plan,
							  int* planlength)
{
	// Generate sampling tool
	std::uniform_real_distribution<double> goal_bias_generation(0.0, 1.0);
	std::uniform_real_distribution<double> joint_ang_generation(0.0, 2*PI);
	std::random_device rd;
	std::default_random_engine generator(rd());

	int epsilon = 30; //epsilon here specifies max. num of steps for each iteration
	double step_size = PI/180;

	int num_vertices = 0;

	// Create tree objects
	RRTTree tree_s(numofDOFs);
	RRTTree tree_g(numofDOFs);

	// Add initial pose as vertex 
	std::vector<double> start_config(armstart_anglesV_rad, armstart_anglesV_rad+numofDOFs);
	std::vector<double> goal_config(armgoal_anglesV_rad, armgoal_anglesV_rad+numofDOFs);
	tree_s.addVertex(start_config);
	tree_g.addVertex(goal_config);
	num_vertices += 2;

	// Initialize new / sampled configurations
	std::vector<double> new_config(numofDOFs);
	std::vector<double> rand_config(numofDOFs);
	// Initialize vector to store nearest neighbor's configuration
	std::vector<double> nn_config(numofDOFs);
	int nn_index = -1;

	// Create while loop until the two trees are connected
	bool connected = false;
	int num_iterations = 0; 
	
	while (!connected) {
		// if (num_iterations % 2 == 0)->Extend tree_s to q_rand; Extend tree_g to q_new 
		// if (num_iterations % 2 == 1)->Extend tree_g to q_rand; Extend tree_s to q_new 

		// Introduce goal bias here
		double bias = goal_bias_generation(generator);
		if (bias > 0.90) {	
			// Set rand_config to goal configuration
			if (num_iterations % 2 == 0)
				copy(armgoal_anglesV_rad, armgoal_anglesV_rad+numofDOFs, rand_config.begin());
			else 
				copy(armstart_anglesV_rad, armstart_anglesV_rad+numofDOFs, rand_config.begin());
		} 
		else {
			// Sample new arm configuration 
			for (int i = 0; i < numofDOFs; i++)
			rand_config[i] = joint_ang_generation(generator);
		}
		// Find nearest neighbor 
		if (num_iterations % 2 == 0) {
			nn_index = tree_s.getNearestVertex(rand_config);
			nn_config = tree_s.getNodeConfig(nn_index);
		} else {
			nn_index = tree_g.getNearestVertex(rand_config);
			nn_config = tree_g.getNodeConfig(nn_index);
		} 

		// Ignore sample that is too close to existing vertices
		if (bias < 0.90 && tree_s.calculateDistance(rand_config, nn_config) < 0.20) continue;

		// Attempt to extend to new sample 
		bool extend = true;
		double* temp_config = &nn_config[0];
		double* backup_config = (double*) malloc(numofDOFs*sizeof(double));
		int step_count = 0;

		while (extend) {
			// save configuration before increment 
			copy(temp_config, temp_config+numofDOFs, backup_config);
			// increment on each joint angle
			incrementToJointAngles(numofDOFs, step_size, &temp_config, rand_config);

			// check collision 
			if (!IsValidArmConfiguration(temp_config, numofDOFs, map, x_size, y_size)) {
				copy(backup_config, backup_config+numofDOFs, temp_config);
				break;
			}
			step_count++;
			if (step_count > epsilon) extend = false;
		}

		// If not trapped...
		if (step_count > 0) {
			// Update ACUTUAL new configuration
			copy(temp_config, temp_config+numofDOFs, new_config.begin());
			
			if (num_iterations % 2 == 0) {
				// Add new_config to the tree
				tree_s.addVertex(new_config);
				tree_s.addEdge(nn_index, tree_s.getNodeID()-1);

				// Find nearest neighbor of new_config from tree_g
				nn_index = tree_g.getNearestVertex(new_config);
				nn_config = tree_g.getNodeConfig(nn_index);

			} else {
				// Add new_config to the tree
				tree_g.addVertex(new_config);
				tree_g.addEdge(nn_index, tree_g.getNodeID()-1);
				// Find nearest neighbor of new_config from tree_s
				nn_index = tree_s.getNearestVertex(new_config);
				nn_config = tree_s.getNodeConfig(nn_index);
			}
			num_vertices++;
			if (num_vertices % 100 == 0) std::cout << "Num of Vertices: " << num_vertices << std::endl;

			// Attempt to extend to new sample 
			// Reset parameters
			extend = true;
			temp_config = &nn_config[0];
			step_count = 0;

			while (extend) {
				// save configuration before increment 
				copy(temp_config, temp_config+numofDOFs, backup_config);
				// increment on each joint angle
				incrementToJointAngles(numofDOFs, step_size, &temp_config, new_config);

				step_count++;
				if (step_count > epsilon) extend = false;
				
				// check collision 
				if (!IsValidArmConfiguration(temp_config, numofDOFs, map, x_size, y_size)) {
					copy(backup_config, backup_config+numofDOFs, temp_config);
					extend = false;
				}
				// check if two trees are connected
				double dist_to_new = 0.0;
				for (int i = 0; i < numofDOFs; i++) {
					double joint_dist = fabs(new_config[i] - temp_config[i]);
					dist_to_new += (joint_dist > PI) ? (2*PI-joint_dist) : joint_dist;
				}
				if (dist_to_new < 0.1) {
					connected = true;
					extend = false;
				} 
			}

			if (step_count > 0) {
				// Update ACUTUAL new configuration
				if (!connected) copy(temp_config, temp_config+numofDOFs, new_config.begin());

				// Add new_config to the tree
				if (num_iterations % 2 == 0) {
					tree_g.addVertex(new_config);
					tree_g.addEdge(nn_index, tree_g.getNodeID()-1);					
				} else {
					tree_s.addVertex(new_config);
					tree_s.addEdge(nn_index, tree_s.getNodeID()-1);
				}
				num_vertices++;
				if (num_vertices % 100 == 0) std::cout << "Num of Vertices: " << num_vertices << std::endl;
			}
		}
		// remember to free backup pointer
		free(backup_config);

		num_iterations++;
	}

	std::cout << "generating the plan" << std::endl; 
	//back-track the path and return the plan 
	std::vector<int> plan1_ids;
	std::vector<int> plan2_ids;
	if (num_iterations % 2 == 0) {
		plan1_ids = tree_s.returnPlan(nn_index);
		plan2_ids = tree_g.returnPlan();
	} else {
		plan1_ids = tree_s.returnPlan();
		plan2_ids = tree_g.returnPlan(nn_index);	
	}
	int path_length = plan1_ids.size() + plan2_ids.size();

	double total_cost = 0.0;
	std::vector<double> config;
	std::vector<double> prev_config;
	*plan = (double**) malloc(path_length*sizeof(double*));
	for (int i = 0; i < plan1_ids.size(); i++) {
		(*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
		config = tree_s.getNodeConfig(plan1_ids[plan1_ids.size()-i-1]);
		copy(config.begin(), config.end(), (*plan)[i]);

		if (i > 0) total_cost += calculateCost(prev_config, config);
		prev_config = config;
	}
	prev_config = tree_g.getNodeConfig(plan1_ids[0]);
	for (int i = 0; i < plan2_ids.size(); i++) {
		(*plan)[i+plan1_ids.size()] = (double*) malloc(numofDOFs*sizeof(double)); 
		config = tree_g.getNodeConfig(plan2_ids[i]);
		copy(config.begin(), config.end(), (*plan)[i+plan1_ids.size()]);

		total_cost += calculateCost(prev_config, config);
		prev_config = config;
	}
	std::cout << "Cost to goal: " << total_cost << std::endl;
	*planlength = path_length;
	
	return;
}

bool isPathValid(std::vector<double> new_config, std::vector<double>& neighbor_config, 
	             int numofDOFs, double* map, int x_size, int y_size)
{
	double distance = 0;
    for (int i = 0; i < numofDOFs; i++){
        if(distance < fabs(new_config[i] - neighbor_config[i]))
            distance = fabs(new_config[i] - neighbor_config[i]);
    }
    int numofsamples = (int)(distance/(PI/90));
    if (numofsamples < 2) return true;

    double* temp_config = &new_config[0];
    for (int i = 0; i < numofsamples; i++) {
    	for (int i = 0; i < numofDOFs; i++) {
    		temp_config[i] += (neighbor_config[i] - new_config[i])/(double)numofsamples;
    	}
    	if (!IsValidArmConfiguration(temp_config, numofDOFs, map, x_size, y_size)) return false;
	}
	return true;
}

static void plannerRRTStar(double*	map, 
						   int x_size,
						   int y_size,
						   double* armstart_anglesV_rad,
						   double* armgoal_anglesV_rad,
						   int numofDOFs,
						   double*** plan,
						   int* planlength)
{
	std::uniform_real_distribution<double> goal_bias_generation(0.0, 1.0);
	std::uniform_real_distribution<double> joint_ang_generation(0.0, 2*PI);
	std::random_device rd;
	std::default_random_engine generator(rd());

	int epsilon = 60; //epsilon here specifies max. num of steps for each iteration
	double step_size = PI/90;

	int num_vertices = 0;

	// Create tree object 
	RRTTree tree(numofDOFs);

	// Add initial pose as vertex 
	std::vector<double> start_config(armstart_anglesV_rad, armstart_anglesV_rad+numofDOFs);
	tree.addVertex(start_config);
	tree.setVertexCost(tree.getNodeID()-1, 0.0);
	num_vertices++;

	// Create while loop until goal configuration is reached 
	bool goal_reached = false;
	
	while (!goal_reached) {
		// Initialize new/ sampled configuration
		std::vector<double> new_config(numofDOFs);

		// Introduce goal bias here
		double bias = goal_bias_generation(generator);
		if (bias > 0.80) {	
			// Set new_config to goal configuration
			copy(armgoal_anglesV_rad, armgoal_anglesV_rad+numofDOFs, new_config.begin());
		} 
		else {
			// Sample new arm configuration 
			for (int i = 0; i < numofDOFs; i++)
				new_config[i] = joint_ang_generation(generator);
		}

		// Find nearest neighbor 
		int nn_index = tree.getNearestVertex(new_config);
		std::vector<double> nn_config = tree.getNodeConfig(nn_index);

		if (bias < 0.80 && tree.calculateDistance(new_config, nn_config) < 0.1) continue;

		// Attempt to extend to new sample 
		bool extend = true;
		double* temp_config = &nn_config[0];
		double* backup_config = (double*) malloc(numofDOFs*sizeof(double));
		int step_count = 0;

		while (extend) {
			// save configuration before increment 
			copy(temp_config, temp_config+numofDOFs, backup_config);
			// increment on each joint angle
			incrementToJointAngles(numofDOFs, step_size, &temp_config, new_config);

			// check collision 
			if (!IsValidArmConfiguration(temp_config, numofDOFs, map, x_size, y_size)) {
				copy(backup_config, backup_config+numofDOFs, temp_config);
				extend = false;
			} 
			else {
				step_count++;
				if (step_count > epsilon) extend = false;
			} 
		}

		free(backup_config);

		if (step_count > 0) {
			// Update ACUTUAL new configuration
			copy(temp_config, temp_config+numofDOFs, new_config.begin());
			// Add new_config to the tree
			tree.addVertex(new_config);
			num_vertices++;
			
			if (num_vertices % 1000 == 0) std::cout << "Num of Vertices: " << num_vertices << std::endl;

			// Add cost of new vertex
			nn_config = tree.getNodeConfig(nn_index);
			double cost_new = tree.getVertexCost(nn_index) + calculateCost(nn_config, new_config);
			tree.setVertexCost(tree.getNodeID()-1, cost_new);
			
			// Find all neighbors within certain distance from new_config 
			double radius = 1.6*pow(log(num_vertices)/num_vertices, 1/numofDOFs); //epsilon * step_size;
			std::vector<int> near_neighbors = tree.getNearVertices(nn_index, radius);

			//if (num_vertices % 1000 == 0) std::cout << "Num of near neighbors: " << near_neighbors.size() << std::endl;
			
			// Try to find best parent for new config 
			// (Is there better path to get ot new vertex from existing neighbors)
			std::vector<double> neighbor_config;
			int min_neighbor_id = nn_index;
			// record indices of valid neighbors of new config (i.e. collision-free)
			std::vector<int> valid_neighbor_id;
			
			for (int i = 0; i < near_neighbors.size(); i++) {
				neighbor_config = tree.getNodeConfig(near_neighbors[i]);
				if (isPathValid(new_config, neighbor_config, numofDOFs, map, x_size, y_size)) {
					double cost_temp = tree.getVertexCost(near_neighbors[i]) + calculateCost(new_config, neighbor_config);
					if (cost_temp + 0.1 < cost_new) {
						// record index of x_min
						min_neighbor_id = near_neighbors[i];
						// update cost of new vertex
						cost_new = cost_temp;
						tree.setVertexCost(tree.getNodeID()-1, cost_new);
					}
					valid_neighbor_id.push_back(near_neighbors[i]);
				} 
			}

			// Add edge for new vertex now 
			tree.addEdge(min_neighbor_id, tree.getNodeID()-1);

			// Try to update cost of each vertex 
			// (Is there better path (through new vertex) to some existing neighbor w/lower cost)
			for (int i = 0; i < valid_neighbor_id.size(); i++) {
				if (valid_neighbor_id[i] != min_neighbor_id) {
					neighbor_config = tree.getNodeConfig(valid_neighbor_id[i]);
					double cost_temp = cost_new + calculateCost(new_config, neighbor_config);
					// In case cost(near) > cost(new) + cost(new,near)
					// If so, disconnect existing node and its parent AND add new edges 
					if (cost_temp + 0.1 < tree.getVertexCost(valid_neighbor_id[i])) {
						// update neighbor's cost
						tree.setVertexCost(valid_neighbor_id[i], cost_temp);
						// remove connection (edge) between the neighbor and its parent
						tree.removeEdge(valid_neighbor_id[i]);
						// add connection (edge) between the neighbor and new vertex
						tree.addEdge(tree.getNodeID()-1, valid_neighbor_id[i]);
					}
				}
			}
			// Check whether we reach the goal 
			double dist_to_goal = 0.0;
			for (int i = 0; i < numofDOFs; i++) {
				double joint_dist = fabs(new_config[i] - armgoal_anglesV_rad[i]);
				dist_to_goal += (joint_dist > PI) ? (2*PI-joint_dist) : joint_dist;
			}
			if (dist_to_goal < 0.1) goal_reached = true; 
		}
		if (num_vertices > 50000) {
			std::cout << "Aborted due to timeout ... Please re-plan" << std::endl;
			break; 
		}
	}

	if (num_vertices <= 50000) {
	    std::cout << "generating the plan" << std::endl;
		//back-track the path and return the plan
		std::vector<int> plan_ids = tree.returnPlan();
		int path_length = plan_ids.size();

		*plan = (double**) malloc(path_length*sizeof(double*));
		for (int i = 0; i < path_length; i++) {
			(*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
			std::vector<double> config = tree.getNodeConfig(plan_ids[path_length-i-1]);
			copy(config.begin(), config.end(), (*plan)[i]);
		}
		std::cout << "Cost to goal: " << tree.getVertexCost(tree.getNodeID()-1) << std::endl;
		*planlength = path_length;
	}

	return;
}

static void plannerPRM(double*	map, 
						int x_size,
						int y_size,
						double* armstart_anglesV_rad,
						double* armgoal_anglesV_rad,
						int numofDOFs,
						double*** plan,
						int* planlength)
{
	/********************Building Roadmap************************/
	std::uniform_real_distribution<double> joint_ang_generation(0.0, 2*PI);
	std::random_device rd;
	std::default_random_engine generator(rd());

	// Create graph object 
	PRMGraph graph(numofDOFs);

	// K specifies number of neighbors we try to connect with new vertex
	int K = 10;
	int num_vertices = 0;

	bool extend = true;
	double* new_config = (double*) malloc(numofDOFs*sizeof(double));

	while (extend) {
		// Random sampling for new vertex/configuration 
		for (int i = 0; i < numofDOFs; i++) new_config[i] = joint_ang_generation(generator);

		// Check if new vertex/configuration is valid (collision-free)
		if (IsValidArmConfiguration(new_config, numofDOFs, map, x_size, y_size)) {
			// Add new vertex to graph
			std::vector<double> new_vertex(5,0.0);
			copy(new_config, new_config+numofDOFs, new_vertex.begin());
			graph.addVertex(new_vertex);
			num_vertices++;

			if (num_vertices > 100) {
				// Find neighbors of new vertex in the graph
				std::vector<int> knn_id = graph.findKNN(new_vertex, K);

				// Check if the path from new vertex to the neighbor is valid (collision-free)
				std::vector<double> neighbor_vertex;
				for (int i = 0; i < knn_id.size(); i++) {
					neighbor_vertex = graph.getNodeConfig(knn_id[i]);
					// If valid, add edge between new vertex and the neighbor in the graph 
					if (isPathValid(new_vertex, neighbor_vertex, numofDOFs, map, x_size, y_size))
						graph.addEdge(knn_id[i], graph.getCurrentNodeID()-1);
				}
			}
		}
		//if (num_vertices%100 == 0) std::cout << num_vertices << std::endl;
		if (num_vertices >= 1000) extend = false;
	}

	free(new_config);

	std::cout << "Roadmap construction is completed. Searching for path now." << std::endl;

	/********************Searching Path*************************/

	std::vector<double> start_config(armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs);
	std::vector<double> goal_config(armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs); 

	int start_on_graph_id = graph.getNearestVertex(start_config);
	int goal_on_graph_id = graph.getNearestVertex(goal_config);

	std::vector<double> start_on_graph = graph.getNodeConfig(start_on_graph_id);
	std::vector<double> goal_on_graph =  graph.getNodeConfig(goal_on_graph_id);

	// Perform Depth-First Search
	std::stack<int> DFS_stack;
	DFS_stack.push(start_on_graph_id);

    bool goal_reached = false;
	if (start_on_graph_id == goal_on_graph_id) goal_reached = true;

	// Create dictionary with key = node_id (visited) and value = parent_id
	std::unordered_map<int, int> node_info;

	int curr_id;
	while (!goal_reached && !DFS_stack.empty()) {
		// Get the vertex at the top of stack
		curr_id = DFS_stack.top();
		DFS_stack.pop();

		vector<int> curr_neighbors = graph.getNeighborsID(curr_id);
		// Add neighbors to the stack
		for (int i = 0; i < curr_neighbors.size(); i++) {
			if (curr_neighbors[i] == goal_on_graph_id) {
				node_info[goal_on_graph_id] = curr_id;
				goal_reached = true;
				break;
			} 
			else if (node_info.find(curr_neighbors[i]) == node_info.end()) {
				node_info[curr_neighbors[i]] = curr_id;
				DFS_stack.push(curr_neighbors[i]);
			}
		}
	}

	if (goal_reached) {
		std::vector<int> plan_ids;
		int index_iter = goal_on_graph_id;
		double total_cost = 0.0;

		while (index_iter != start_on_graph_id) {
			plan_ids.push_back(index_iter);
			index_iter = node_info[index_iter];
		}
		plan_ids.push_back(start_on_graph_id);

		int path_length = plan_ids.size()+2;
		*plan = (double**) malloc(path_length*sizeof(double*));

		(*plan)[0] = (double*) malloc(numofDOFs*sizeof(double));
		copy(armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs, (*plan)[0]);

		std::vector<double> config;
		std::vector<double> prev_config = start_config;

		for (int i = 0; i < plan_ids.size(); i++) {
			(*plan)[i+1] = (double*) malloc(numofDOFs*sizeof(double)); 
			config = graph.getNodeConfig(plan_ids[plan_ids.size()-i-1]);
			copy(config.begin(), config.end(), (*plan)[i+1]);

			if (i > 0) total_cost += calculateCost(prev_config, config);
			prev_config = config;
		}

		(*plan)[path_length-1] = (double*) malloc(numofDOFs*sizeof(double));
		copy(armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs, (*plan)[path_length-1]);

		total_cost += calculateCost(prev_config, goal_config);

		std::cout << "Cost to goal: " << total_cost << std::endl;
		*planlength = path_length;
	}
	else {
		std::cout << "Unable to connect START and GOAL on exisiting roadmap" << std::endl;
	}


	/**************For visualizing vertices only***************/

	// *plan = (double**) malloc(num_vertices*sizeof(double*));

	// for (int i = 0; i < num_vertices; i++) {
	// 	(*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
	// 	std::vector<double> config = graph.getNodeConfig(i);
	// 	copy(config.begin(), config.end(), (*plan)[i]);
	// 	for(int j = 0; j < numofDOFs; j++) std::cout << (*plan)[i][j] << "  ";
	// 	std::cout << std::endl;
	// }
	// *planlength = num_vertices;
	/*********************************************************/

	return;
}

//prhs contains input parameters (3): 
//1st is matrix with all the obstacles
//2nd is a row vector of start angles for the arm 
//3nd is a row vector of goal angles for the arm 
//plhs should contain output parameters (2): 
//1st is a 2D matrix plan when each plan[i][j] is the value of jth angle at the ith step of the plan
//(there are D DoF of the arm (that is, D angles). So, j can take values from 0 to D-1
//2nd is planlength (int)
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[])
     
{ 
    
    /* Check for proper number of arguments */    
    if (nrhs != 4) { 
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Four input arguments required."); 
    } else if (nlhs != 2) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required."); 
    } 
        
    /* get the dimensions of the map and the map matrix itself*/     
    int x_size = (int) mxGetM(MAP_IN);
    int y_size = (int) mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);
    
    /* get the start and goal angles*/     
    int numofDOFs = (int) (MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
    if(numofDOFs <= 1){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "it should be at least 2");         
    }
    double* armstart_anglesV_rad = mxGetPr(ARMSTART_IN);
    if (numofDOFs != MAX(mxGetM(ARMGOAL_IN), mxGetN(ARMGOAL_IN))){
        	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "numofDOFs in startangles is different from goalangles");         
    }
    double* armgoal_anglesV_rad = mxGetPr(ARMGOAL_IN);
 
    //get the planner id
    int planner_id = (int)*mxGetPr(PLANNER_ID_IN);
    if(planner_id < 0 || planner_id > 3){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidplanner_id",
                "planner id should be between 0 and 3 inclusive");         
    }
    
    //call the planner
    double** plan = NULL;
    int planlength = 0;
    
    //you can may be call the corresponding planner function here
    if (planner_id == RRT)
    {
    	clock_t tStart = clock();
    	plannerRRT(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    	std::cout << "Time taken for planning : " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << std::endl;
    }
    else if (planner_id == RRTCONNECT)
    {
    	clock_t tStart = clock();
    	plannerRRTConnect(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    	std::cout << "Time taken for planning : " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << std::endl;
    }
    else if (planner_id == RRTSTAR)
    {
    	clock_t tStart = clock();
    	plannerRRTStar(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    	std::cout << "Time taken for planning : " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << std::endl;
    }
    else if (planner_id == PRM)
    {
    	clock_t tStart = clock();
    	plannerPRM(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    	std::cout << "Time taken for planning : " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << std::endl;
    }
    else {
    	//dummy planner which only computes interpolated path
    	planner(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength); 
    }
    
    printf("planner returned plan of length=%d\n", planlength); 
    
    /* Create return values */
    if(planlength > 0)
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)planlength, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);        
        //copy the values
        int i,j;
        for(i = 0; i < planlength; i++)
        {
            for (j = 0; j < numofDOFs; j++)
            {
                plan_out[j*planlength + i] = plan[i][j];
            }
            free(plan[i]);
        }
        free(plan);
    }
    else
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);
        //copy the values
        int j;
        for(j = 0; j < numofDOFs; j++)
        {
                plan_out[j] = armstart_anglesV_rad[j];
        }     
    }
    PLANLENGTH_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL); 
    double* planlength_out = mxGetPr(PLANLENGTH_OUT);
    //char* planlength_out = (char*)mxGetData(PLANLENGTH_OUT);
    *planlength_out = planlength;

    
    return;
    
}





