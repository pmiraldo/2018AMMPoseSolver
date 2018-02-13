/*
 * @file SURF_FlannMatcher
 * @brief SURF detector + descriptor + FLANN Matcher
 * @author A. Huaman
 */

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include <Eigen/Dense>
#include <opengv/amm.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/NoncentralRelativeAdapter.hpp>
#include <opengv/optimization_tools/objective_function_tools/GlobalPnPFunctionInfo.hpp>
#include <opengv/optimization_tools/objective_function_tools/SquaredFunctionNoIterationsInfo.hpp>
#include <opengv/optimization_tools/solver_tools/SolverToolsNoncentralRelativePose.hpp>

#include "opencv2/core.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/xfeatures2d.hpp"

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"


using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;
//void readme();
/*
 * @function main
 * @brief Main function
 */
void returnCorrespondences(const Mat & img1, const Mat & img2, std::vector<KeyPoint> & k1, std::vector<KeyPoint> & k2);

int main( int argc, char** argv )
{
  if( argc != 5 )
  {
    std::cout << "argc: " << argc << std::endl;
    //readme();
    return -1;
  }
  Mat img_left_1 = imread( argv[1], IMREAD_GRAYSCALE );
  Mat img_left_2 = imread( argv[2], IMREAD_GRAYSCALE );

  Mat img_right_1 = imread( argv[3], IMREAD_GRAYSCALE );
  Mat img_right_2 = imread( argv[4], IMREAD_GRAYSCALE );
  if( !img_left_1.data || !img_left_2.data || !img_right_1.data || !img_right_2.data)
  {
    std::cout<< " --(!) Error reading images " << std::endl;
    return -1;
  }
  std::vector<KeyPoint> k1_left;
  std::vector<KeyPoint> k2_left;
  returnCorrespondences(img_left_1, img_left_2, k1_left, k2_left);

  std::vector<KeyPoint> k1_right;
  std::vector<KeyPoint> k2_right;
  returnCorrespondences(img_right_1, img_right_2, k1_right, k2_right);

  /*std::cout << "Left: " << std::endl;
  for(int i=0; i < k1_left.size(); ++i){
    std::cout << "(" << k1_left[i].pt.x << ", " << k1_left[i].pt.y << ")  (" << k2_left[i].pt.x << ", " << k2_left[i].pt.y << ")" <<  std::endl;
  }

  std::cout << "Right: " << std::endl;
  for(int i=0; i < k1_right.size(); ++i){
    std::cout << "(" << k1_right[i].pt.x << ", " << k1_right[i].pt.y << ")  (" << k2_right[i].pt.x << ", " << k2_right[i].pt.y << ")" <<  std::endl;
    }*/

   
  std::cout << "Cal matrices: " << std::endl;
  Eigen::Matrix3d K_left = Eigen::Matrix3d::Identity(3,3);
  K_left(0,0) = 9.597910e+02; K_left(0,1) =  0.000000e+00; K_left(0,2) = 6.960217e+02;
  K_left(1,0) = 0.000000e+00; K_left(1,1) =  9.569251e+02; K_left(1,2) = 2.241806e+02;
  K_left(2,0) = 0.000000e+00; K_left(2,1) =  0.000000e+00; K_left(2,2) = 1.000000e+00;

  Eigen::Matrix3d K_right = Eigen::Matrix3d::Identity(3,3);
  K_right(0,0) = 9.037596e+02; K_right(0,1) =  0.000000e+00; K_right(0,2) = 6.957519e+02;
  K_right(1,0) = 0.000000e+00; K_right(1,1) =  9.019653e+02; K_right(1,2) = 2.242509e+02;
  K_right(2,0) = 0.000000e+00; K_right(2,1) =  0.000000e+00; K_right(2,2) = 1.000000e+00;

  //We have the calibration matrices and the vector points now it is necessary to build the bearing vectors and convert them in the reference frame of the centra

  //Pose for viewpoint left
  opengv::translation_t position_left = Eigen::MatrixXd::Zero(3,1);
  position_left(0,0) = -0.27;
  opengv::rotation_t rotation_left    = Eigen::MatrixXd::Identity(3,3);

  //Pose for viewpoint right
  opengv::translation_t position_right = Eigen::MatrixXd::Zero(3,1);
  position_right(0,0) = 0.27;
  opengv::rotation_t rotation_right    = Eigen::MatrixXd::Identity(3,3);

  //Create a fake central camera
 
  opengv::bearingVectors_t bearingVectors1;
  opengv::bearingVectors_t bearingVectors2;
  Eigen::MatrixXd aux = Eigen::MatrixXd::Zero(3,1);
  for(int j = 0; j < k1_left.size(); ++j){
    aux(0,0) = k1_left[j].pt.x; aux(1,0) = k1_left[j].pt.y; aux(2,0) = 1;
    aux = K_left.inverse() * aux;
    bearingVectors1.push_back(aux);
  }
  
  for(int j = 0; j < k1_right.size(); ++j){
    aux(0,0) = k1_right[j].pt.x; aux(1,0) = k1_right[j].pt.y; aux(2,0) = 1;
    aux = K_right.inverse() * aux;
    bearingVectors1.push_back(aux);
  }
  
  for(int j = 0; j < k2_left.size(); ++j){
    aux(0,0) = k2_left[j].pt.x; aux(1,0) = k2_left[j].pt.y; aux(2,0) = 1;
    aux = K_left.inverse() * aux;
    bearingVectors2.push_back(aux);
  }
  
  for(int j = 0; j < k2_right.size(); ++j){
    aux(0,0) = k2_right[j].pt.x; aux(1,0) = k2_right[j].pt.y; aux(2,0) = 1;
    aux = K_right.inverse() * aux;
    bearingVectors2.push_back(aux);
  }
  std::cout << "bearing vector1: " << bearingVectors1.size() << std::endl;
  std::cout << "bearing vector2: " << bearingVectors2.size() << std::endl;
  //To use the methods it will be necessary to build the central adapter
  /*opengv::bearingVectors_t v1;
  opengv::bearingVectors_t v2;
  int cols = (int)good_matches.size();
  for( int i = 0; i < cols; i++ )
    {
      std::cout << "Start cycle: " << i << std::endl;
      Eigen::MatrixXd v_left = Eigen::MatrixXd::Zero(3,1);
      v_left(0,0) = k1[i].pt.x;
      v_left(1,0) = k1[i].pt.y;
      v_left(2,0) = 1;
      //std::cout << "v1: " << std::endl << v1 << std::endl;

      Eigen::MatrixXd v_right = Eigen::MatrixXd::Zero(3,1);
      v_right(0,0) = k2[i].pt.x;
      v_right(1,0) = k2[i].pt.y;
      v_right(2,0) = 1;
      //std::cout << "v2: " << std::endl << v2 <<std::endl;

      v_left = K_left.inverse() * v_left;
      double norm_left = v_left.norm();
      v_left = v_left / norm_left;

      v_right = K_right.inverse() * v_right;
      double norm_right = v_right.norm();
      v_right = v_right /norm_right;
      
      v1.push_back(v_left);
      v2.push_back(v_right);
    }
  opengv::rotations_t camRotations;
  opengv::rotation_t rotation1 = Eigen::MatrixXd::Identity(3,3);
  opengv::rotation_t rotation2 = Eigen::MatrixXd::Identity(3,3);
  opengv::translation_t position1 = Eigen::MatrixXd::Zero(3,1);
  opengv::translation_t position2 = Eigen::MatrixXd::Zero(3,1);
  position2(0,0) = 0.54;
  camRotations.push_back(rotation1);
  camRotations.push_back(rotation2);
 
  opengv::translations_t camOffsets;
  camOffsets.push_back(position1);
  camOffsets.push_back(position2);
  std::vector<int> camCorrespondences1;
  std::vector<int> camCorrespondences2;
  int camCorrespondence = 0;
  for(int i = 0; i < good_matches.size(); ++i){
    camCorrespondences1.push_back(camCorrespondence);
    camCorrespondences2.push_back(camCorrespondence++);
  }
  
  //Extract the relative pose
  opengv::translation_t position; opengv::rotation_t rotation;
  opengv::extractRelativePose(position1, position2, rotation1, rotation2, position, rotation, false );

  //create non-central relative adapter
  opengv::relative_pose::NoncentralRelativeAdapter adapter(
						       v1,
						       v2,
						       camCorrespondences1,
						       camCorrespondences2,
						       camOffsets,
						       camRotations,
						       position,
						       rotation);
  return 0;*/
}

/*
 * @function readme
 */
/*void readme(){
  std::cout << " Usage: ./SURF_FlannMatcher <img1> <img2>" << std::endl;
  }*/

void returnCorrespondences(const Mat & img1, const Mat & img2, std::vector<KeyPoint> & k1, std::vector<KeyPoint> & k2){
  //-- Step 1: Detect the keypoints using SURF Detector, compute the descriptors
  int minHessian = 400;
  Ptr<SURF> detector = SURF::create();
  detector->setHessianThreshold(minHessian);
  
  std::vector<KeyPoint> keypoints_1, keypoints_2;
  Mat descriptors_1, descriptors_2;
  detector->detectAndCompute( img1, Mat(), keypoints_1, descriptors_1 );
  detector->detectAndCompute( img2, Mat(), keypoints_2, descriptors_2 );
  //-- Step 2: Matching descriptor vectors using FLANN matcher
  FlannBasedMatcher matcher;
  std::vector< DMatch > matches;
  matcher.match( descriptors_1, descriptors_2, matches );
  double max_dist = 0; double min_dist = 100;
  //-- Quick calculation of max and min distances between keypoints
  for( int i = 0; i < descriptors_1.rows; i++ )
  {
    double dist = matches[i].distance;
    if( dist < min_dist ) min_dist = dist;
    if( dist > max_dist ) max_dist = dist;
  }
  printf("-- Max dist : %f \n", max_dist );
  printf("-- Min dist : %f \n", min_dist );
  //-- Draw only "good" matches (i.e. whose distance is less than 2*min_dist,
  //-- or a small arbitary value ( 0.02 ) in the event that min_dist is very
  //-- small)
  //-- PS.- radiusMatch can also be used here.
  std::vector< DMatch > good_matches;
  for( int i = 0; i < descriptors_1.rows; i++ )
  {
    if( matches[i].distance <= max(2*min_dist, 0.02) ){
      good_matches.push_back( matches[i]);
    }
  }
  //-- Draw only "good" matches
  k1.clear();
  k2.clear();
  for( int i = 0; i < (int)good_matches.size(); i++ )
  {
    //printf( "-- Good Match [%d] Keypoint 1: %d  -- Keypoint 2: %d  \n", i, good_matches[i].queryIdx, good_matches[i].trainIdx );
    std::cout << ( "-- Good Match [%d] Keypoint 1: %d  -- Keypoint 2: %d  \n", i, good_matches[i].queryIdx, good_matches[i].trainIdx ) << std::endl;
    k1.push_back(keypoints_1[good_matches[i].queryIdx]);
    k2.push_back(keypoints_2[good_matches[i].trainIdx]);    
  }
  
}
