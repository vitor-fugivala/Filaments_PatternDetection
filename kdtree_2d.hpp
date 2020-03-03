//=================================
//       Include Guard           //
//=================================
#ifndef KDTREE2_HPP
#define KDTREE2_HPP

//=================================
//       Included Dependencies   //
//=================================
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <algorithm>

//=================================
//       Struct Definations      //
//=================================
using namespace std;

struct point {
    double x,y;
    int subset_no; //for the jacknife
    double radius;
};

struct box {
    double min_x,max_x,min_y,max_y;
};

//contains a list of the points sorted by x,y,z
//invariant: x.size() == y.size() == z.size()
//invariant: is_sorted(x), is_sorted(y), is_sorted(z)
struct sorted_points {
    vector<point> x,y;
};

struct kd_node {
    vector<point> points; //only Leafs will store points
    box bounding_box;
    kd_node *right;
    kd_node *left;   
};

typedef struct kd_node* kd_tree;

//=================================
//       KD-Tree Functions       //
//=================================

//Requires: pts has at least one point
//Ensures: return a kd-tree of the points, 
//         cuts in plane of largest dimension
kd_tree new_kdtree(vector<point> points);

//Free Memory of kdtree
void free_kdtree(kd_tree T);

//Prints nodes and leaves
void print_tree(kd_tree T, int depth);

//Returns true if T is a leaf node and contains points
bool is_leaf(kd_tree T);

//Computes minumum distance between two nodes
//Requires: T1 != NULL, T2 != NULL
double min_distance(kd_tree T1, kd_tree T2);

//Computes maximum distance between two nodes
//Requires: T1 != NULL, T2 != NULL
double max_distance(kd_tree T1, kd_tree T2);

#endif






