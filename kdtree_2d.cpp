//=================================
//     KD-Tree Implementation    //
//=================================
#include "kdtree_2d.hpp"
using namespace std;

//Smallest number of data points in a leaf
#define MIN_ELEMS 5 

//returns largest dimension of sorted_points structure
int largest_dim(sorted_points &P){
    double dx = P.x.back().x - P.x.front().x;
    double dy = P.y.back().y - P.y.front().y;;
    if (dx >= dy ) {return 0;}
    else return 1;
}

//splits the points by the median of the largest dimension
//requires that there are at least 2 points
vector<sorted_points> split(sorted_points &P){
    vector<sorted_points> split_vectors;
    int dim = largest_dim(P);
    sorted_points left,right;

    vector<point> xpts = P.x;
    vector<point> ypts = P.y;

    if (dim == 0){
        size_t m = xpts.size()/2;
        double x_cut = xpts[m].x;
        for (size_t i = 0; i < xpts.size(); i++){
            if (i < m) left.x.push_back(xpts[i]);
            else right.x.push_back(xpts[i]);
          
            if (ypts[i].x < x_cut) left.y.push_back(ypts[i]);
            else right.y.push_back(ypts[i]);
        }
    }
    
    else {
        size_t m = ypts.size()/2;
        double y_cut = ypts[m].y;
        for (size_t i = 0; i < ypts.size(); i++){
            if (i < m) left.y.push_back(ypts[i]);
            else right.y.push_back(ypts[i]);
          
            if (xpts[i].y < y_cut) left.x.push_back(xpts[i]);
            else right.x.push_back(xpts[i]);

        }
    }
    
    split_vectors.push_back(left);
    split_vectors.push_back(right);
    return split_vectors;
}


box find_bounding_box(vector<point> &P){
    box B;
    B.min_x =  99999.9;
    B.max_x = -99999.9;
    B.min_y =  99999.9;
    B.max_y = -99999.9;
    for (size_t i = 0; i < P.size(); i++){
        point pt = P[i];
        if (pt.x < B.min_x) {B.min_x = pt.x;}
        if (pt.x > B.max_x) {B.max_x = pt.x;}
        if (pt.y < B.min_y) {B.min_y = pt.y;}
        if (pt.y > B.max_y) {B.max_y = pt.y;}
    }
    return B;
}

box merge_box(box B1, box B2){
    box B0;
    B0.min_x = min(B1.min_x,B2.min_x);
    B0.max_x = max(B1.max_x,B2.max_x);
    B0.min_y = min(B1.min_y,B2.min_y);
    B0.max_y = max(B1.max_y,B2.max_y);
    return B0;
}

bool cmp_pts_x(point pt1, point pt2){
    return pt1.x < pt2.x;
}
bool cmp_pts_y(point pt1, point pt2){
    return pt1.y < pt2.y;
}



//requires: pts has at least one point, depth >= 0
//ensures: return a kd-tree of the points
kd_node* make_kdtree(sorted_points &P){
    kd_node *T = new kd_node;
    
    // If there are less then MIN_ELEMS points, this is a leaf
    // The factor of 2 is because we can't split otherwise
    if (P.x.size() <  2*MIN_ELEMS) {
        T->points = P.x;
        T->bounding_box = find_bounding_box(P.x);
        T->left = NULL;
        T->right = NULL;
    }
    else{ //Else it is a node and we recursively build the tree
        vector<sorted_points> split_vectors = split(P);
        T->left = make_kdtree(split_vectors[0]);
        T->right = make_kdtree(split_vectors[1]);
        T->bounding_box = merge_box(T->left->bounding_box, 
                                    T->right->bounding_box);
    }
    return T;
}

//Free Memory of kdtree
void free_kdtree(kd_tree T){
    if (T!= NULL && T->left != NULL){
        free_kdtree(T->left);
        free_kdtree(T->right);
    }
    delete(T);
}

//requires: pts has at least one point
//ensures: return a kd-tree of the points, 
//         rotates cuts in x,y,z in that order
kd_tree new_kdtree(vector<point> points){
    sorted_points P; 
    vector<point> x_points (points);
    vector<point> y_points (points);
    sort(x_points.begin(),x_points.end(),cmp_pts_x);
    sort(y_points.begin(),y_points.end(),cmp_pts_y);
    P.x = x_points;
    P.y = y_points;
    return make_kdtree(P);  
}

//Prints nodes and leaves
void print_tree(kd_tree T, int depth){
    if (T->left == NULL){
        cout << "LEAF. Depth: " << depth << endl;
        for (size_t i = 0; i < T->points.size(); i++){
            point pt = T->points[i];
            cout << "point: " << pt.x << "," << pt.y  << endl;
        }
    }
    else{
        cout << "NODE. Depth: " << depth << endl;
        print_tree(T->left, depth + 1);
        print_tree(T->right, depth + 1);
    }
}

//Returns true if T is a leaf node and contains points
bool is_leaf(kd_tree T){
    if (T != NULL && T->left == NULL && T->right == NULL){
        return true;
    }
    return false;
}

//Computes minumum distance between two nodes
//Requires: T1 != NULL, T2 != NULL
double min_distance(kd_tree T1, kd_tree T2){
    double distance = 0.0;
    if (T2->bounding_box.max_x < T1->bounding_box.min_x){
        double diff = T2->bounding_box.max_x - T1->bounding_box.min_x;
        distance += diff * diff;
    }
    else if (T2->bounding_box.min_x > T1->bounding_box.max_x){
        double diff = T2->bounding_box.min_x - T1->bounding_box.max_x;
        distance += diff * diff;
    }
    if (T2->bounding_box.max_y < T1->bounding_box.min_y){
        double diff = T2->bounding_box.max_y - T1->bounding_box.min_y;
        distance += diff * diff;
    }
    else if (T2->bounding_box.min_y > T1->bounding_box.max_y){
        double diff = T2->bounding_box.min_y - T1->bounding_box.max_y;
        distance += diff * diff;
    }
   return sqrt(distance);
}

//Helper function to compute max of 4 doubles
double maxx(double w, double x, double y, double z){
    return max(w,max(x,max(y,z)));
}

//Computes maximum distance between two nodes
//Requires: T1 != NULL, T2 != NULL
double max_distance(kd_tree T1, kd_tree T2){
    double distance = 0.0;
    double diffx = maxx(fabs(T1->bounding_box.min_x - T2->bounding_box.min_x),
                        fabs(T1->bounding_box.min_x - T2->bounding_box.max_x),
                        fabs(T1->bounding_box.max_x - T2->bounding_box.min_x),
                        fabs(T1->bounding_box.max_x - T2->bounding_box.max_x));
    double diffy = maxx(fabs(T1->bounding_box.min_y - T2->bounding_box.min_y),
                        fabs(T1->bounding_box.min_y - T2->bounding_box.max_y),
                        fabs(T1->bounding_box.max_y - T2->bounding_box.min_y),
                        fabs(T1->bounding_box.max_y - T2->bounding_box.max_y));

    distance += diffx*diffx + diffy*diffy ;
    return sqrt(distance);
}
