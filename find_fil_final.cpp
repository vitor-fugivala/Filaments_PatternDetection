//=====================================================================//
// Filament Finding Using SCMS (Subspace Constrained Mean Shift)       //
// --------------------------------------------------------------------//
// Cosmic Web Reconstruction through Density Ridges                    //
// http://arxiv.org/abs/1501.05303v3                                   //
//=====================================================================//
#include <cstdlib>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <Eigen>
#include <getopt.h>
#include "./kdtree_2d.hpp"
#include <string>
using namespace Eigen;
using namespace std;

//=====================================================================//
//                         Helper Functions                            //
//=====================================================================//
//combine 2 vectors
//This might be a slow way of doing this
vector<Vector2d> combine_vectors(vector<Vector2d> &A, vector<Vector2d> &B){
    vector<Vector2d> AB;
    AB.reserve( A.size() + B.size() ); // preallocate memory
    AB.insert( AB.end(), A.begin(), A.end() );
    AB.insert( AB.end(), B.begin(), B.end() );
    return AB;
}

//wrapper for min distance
double min_point_distance(Vector2d pt_vector, kd_tree pt_tree){
    vector<point> pts;
    point pt;
    pt.x = pt_vector(0);
    pt.y = pt_vector(1);

    pts.push_back(pt);
    kd_tree leaf_node = new_kdtree(pts);
    return min_distance(leaf_node,pt_tree);
}


//returns all points in the kd tree whithin 3sigma of supplied point
vector<Vector2d> nearestNeighbours(Vector2d src_point, kd_tree point_tree, 
                                   float sigma){
    vector<Vector2d> result;
    if (is_leaf(point_tree)){
        vector<point> pts = point_tree->points;
        for (size_t i=0; i<pts.size(); i++){
            point pt = pts[i];
            Vector2d candidate_pt(pt.x,pt.y);
            if (((src_point - candidate_pt).norm())<3*sigma){
                result.push_back(candidate_pt);
            }
        }
        return result;
    }
    else if (min_point_distance(src_point,point_tree) > 3*sigma){
        return result; //Prune points that are too far away
    }
    else {
        //Recursive case
        vector<Vector2d> left;
        vector<Vector2d> right;
         
        left = nearestNeighbours(src_point, point_tree->left, sigma);
        right = nearestNeighbours(src_point, point_tree->right, sigma);
        
        return combine_vectors(left,right);
    }
}

//Convert data file to vector<Vector2d> and vector<point>
void read_points(string in, vector<Vector2d> *points, 
                 vector<point> *kd_points){
    ifstream infile;
    infile.open(in.c_str());
    string line;
    while(getline(infile,line) && line[0] != '#'){
        Vector2d pt;
        point kd_point;
        stringstream ss(line);
        double x,y;
        ss >> x >> y ;
        pt(0) = x;
        pt(1) = y;
        kd_point.x = x;
        kd_point.y = y;
        points->push_back(pt);
        kd_points->push_back(kd_point);
    }
    infile.close();
}

//Gaussian Kernal value
float std_normal(float x){
    return exp(-0.5*x*x)/sqrt(2*M_PI);
}

//Computes p(x)
float density_estimator(float h, Vector2d x, vector<Vector2d> &points, 
                        kd_tree points_kdtree){ 
    float sum = 0.0;
    //Find nearest nieghbours
    vector<Vector2d> nNs = nearestNeighbours(x,points_kdtree,h);
    for (size_t i = 0; i < nNs.size(); i++){
        sum += std_normal((x-nNs[i]).norm()/h);
    }
    return sum/(points.size()*h*h);
}

//Output a file of thresholded point
void write_output(vector<Vector2d> &points, string suffix){
    stringstream ss;
    ss << suffix;
    cout << "Writing points to: " << ss.str() << endl;
    ofstream fout;
    fout.open(ss.str().c_str(), ios::out);
    for (size_t i = 0; i < points.size(); i++){
        fout << points[i](0) << ',' << points[i](1) << endl;
    }
}

//returns the threshold value as the standard deviation of densities
float calculate_threshold(vector<float> &v){
    float sum = accumulate(v.begin(), v.end(), 0.0);
    float mean = sum / v.size();
    
    vector<float> diff(v.size());
    transform(v.begin(),v.end(),diff.begin(),bind2nd(minus<float>(),mean));
    float sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    float stdev = sqrt(sq_sum / v.size());
    
    return stdev;
}

//=====================================================================//
//                         Main Function                               //
//=====================================================================//
int main(int argc, char* argv[]){
    //Parse Input From Command Line
    cout << "Reading Inputs..." << endl;
    int opt,dflag=0,hflag=0,iflag=0,max_iterations,tflag=0,rflag=0,eflag=0;
    string data_filename,out_filename;
    float h,t,epsilon;
    while ((opt = getopt(argc, argv, "d:h:i:t:e:r:n:")) != -1){
        switch(opt){
            case 'd': // Name of the data file of points
                data_filename = optarg;
                cout << "    Data filename: " << data_filename << endl;
                dflag = 1;
                break;
            case 'h': // Smoothing bandwidth h
                h = float(atof(optarg));
                hflag = 1;
                cout << "    Value of h: " << h << endl;
                break;
            case 'i': // iterations
                max_iterations = atoi(optarg);
                iflag = 1;
                cout << "    Max Iterations: " << max_iterations << endl;
                break;
            case 't': // threshold
                t = atoi(optarg);
                tflag = 1;
                break;
            case 'r': // resolution 
                rflag = 1;
                break;
            case 'e': // epsilon
                epsilon = atof(optarg);
                eflag = 1;
                cout << "    Value of e: " << epsilon << endl;
                break;
            case 'n': // name of output file 
                out_filename = optarg;
                cout << "    Output filename: " << out_filename << endl;
                break;
            case '?':
                cout << "Error in entered arguments" << endl;
                return 1;
            default:
                abort();            
        }
    }
    cout << "Done." << endl;
    
    //Check that all arguments were received
    if ((dflag == 0) || (hflag == 0) || (iflag == 0)){
        cout << "You missed an argument!" << endl;
        return 1;
    }
    
    //Extracting Information from data file - Input Data : {X1,...,Xn}
    cout << "Extracting Information from data file..." ;
    vector<Vector2d> points;
    vector<point> kd_points;
    read_points(data_filename,&points,&kd_points);
    cout << "Done." << endl;
    cout << "Building kd tree of points...";
    kd_tree points_kdtree = new_kdtree(kd_points);
    cout << "Done." << endl;
    cout << "No. of points: " << points.size() << endl;
        
    //Step 1: Thresholding: remove x in X if p(x) < t
    //        calculate t using the RMS in the density fluctuations
    cout << "Thresholding...";
    //Find all the densities
    vector<float> densities(points.size());
    for (size_t i = 0; i < points.size(); i++){
        densities[i] = density_estimator(h,points[i],points,points_kdtree);
    }
    //calculate the threshold value : standard deviation of densities
    if (tflag == 0){
        t = calculate_threshold(densities);
    }
    cout << "The Threshold value is " << t << endl;
    //Then perform thresholding
    vector<Vector2d> thresholded_points;
    thresholded_points.reserve(points.size());
    for (size_t i = 0; i < points.size(); i++){
        if (densities[i] > t){
            thresholded_points.push_back(points[i]);
        }
    }
    cout << "No. of thresholded points: " << thresholded_points.size() << endl;
    
    if (rflag == 1){
        //Optional increase of resolution
        size_t holder = thresholded_points.size();
        for (size_t i = 0; i < holder; i++){
            Vector2d point1 = thresholded_points[i];
            point1(0) += h/2;
            point1(1) += h/2;
            thresholded_points.push_back(point1);
            
            Vector2d point2 = thresholded_points[i];
            point2(0) += h/2;
            point2(1) += h/2;
            thresholded_points.push_back(point2);
            
            Vector2d point3 = thresholded_points[i];
            point3(0) += h/2;
            point3(1) -= h/2;
            thresholded_points.push_back(point3);
            
            Vector2d point4 = thresholded_points[i];
            point4(0) += h/2;
            point4(1) -= h/2;
  
            thresholded_points.push_back(point4);
            
            Vector2d point5 = thresholded_points[i];
            point5(0) -= h/2;
            point5(1) += h/2;
            thresholded_points.push_back(point5);
            
            Vector2d point6 = thresholded_points[i];
            point6(0) -= h/2;
            point6(1) += h/2;
            thresholded_points.push_back(point6);
            
            Vector2d point7 = thresholded_points[i];
            point7(0) -= h/2;
            point7(1) -= h/2;
            thresholded_points.push_back(point7);
            
            Vector2d point8 = thresholded_points[i];
            point8(0) -= h/2;
            point8(1) -= h/2;
            thresholded_points.push_back(point8);
        }
    }    
             
    //Step 2: Peform the following SCMS until convergence
    cout << "Performing SCMS..." << endl;
    if (eflag == 0){
        epsilon = 1.0;
    }
    float max_error = 1.0e14;
    int   iteration = 0;

    //initialize all errors to some arbitrary high value
    ArrayXf errors(thresholded_points.size());
    for (size_t i = 0; i < thresholded_points.size(); i++){
        errors(i) = max_error;
    } 
    
    while ((max_error > epsilon) && (iteration < max_iterations)){
        cout << "    Iteration: " << iteration;
        //For each iteration, we find the maximum error
        max_error = 0.0;
        int points_moved = 0;
        
        for (size_t i = 0; i < thresholded_points.size(); i++){
            //Do nothing for converged pts
            if (errors(i) < epsilon){continue;}     
            points_moved += 1;     
            Vector2d x = thresholded_points[i];
            
            //Compute The Hessian via mu[i] and c[i]
            Matrix2d H = Matrix2d::Zero();
            Vector2d mu(0,0);
            float c = 0.0;          
            Vector2d sum_cx(0,0); //for msv
            float sum_c = 0.0;      //for msv
            
            //Find points that are close enough (whithin 3*h <=> 3 sigma)
            vector<Vector2d> nNs = nearestNeighbours(x,points_kdtree,h);
            if (nNs.size() == 0){
                errors(i) = 0;           
                cout << "Error: Isolated point - should not occur" << endl;
                continue;
            }
            
            for (size_t j = 0; j < nNs.size(); j++){
                mu = (x - nNs[j])/(h*h);
                c  = std_normal((x-nNs[j]).norm()/h);
                H += ((c/points.size())
                     *((mu*mu.transpose())-(MatrixXd::Identity(2,2)/(h*h))));
                
                //for msv
                sum_cx = sum_cx + c*nNs[j];
                sum_c  = sum_c + c;
            }
            
            //Find smallest d-1 eigenvectors
            //Since H is always symmetric, we use an optimized solver
            //with the side benifit the eigenvalues and vectors are sorted 
            //in increasing order
            SelfAdjointEigenSolver<MatrixXd> es(H);
            Matrix2d eigenvectors = es.eigenvectors();
            
            //sort in ascending order and drops the first element/col
            Matrix<double,2,1> ev_red;
            ev_red = (eigenvectors.block<2,1>(0,0)).rowwise().reverse();
            
            //Move the point, msv = mean shift vector
            Vector2d msv = (sum_cx/sum_c) - x;
            Vector2d shift = ev_red*ev_red.transpose()*msv;   
            thresholded_points[i] = thresholded_points[i] + shift;
            
            //Update errors
            float difference = shift.norm();
            errors(i) = difference;
            
            if (difference > max_error){
                max_error = difference;
            }
        }
        cout << "    Points moved: "  << points_moved << endl;
        iteration += 1;
    }
    
    if (iteration == max_iterations) {
        cout << "    Max Iterations Reached." << endl;
    }
    
    cout << "Done." << endl;
    
    //Output: the collection of all remanining points
    write_output(thresholded_points, out_filename);
    free_kdtree(points_kdtree);
    return 0;
}
