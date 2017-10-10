function[armplan] = armplanner(envmap, armstart, armgoal);


 %call the planner in C
 [armplan, armplanlength] = planner(envmap, armstart, armgoal);

