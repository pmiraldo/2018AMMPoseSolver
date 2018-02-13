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
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <opengv/optimization_tools/objective_function_tools/GlobalPnPFunctionInfo.hpp>
#include <opengv/optimization_tools/objective_function_tools/SquaredFunctionNoIterationsInfo.hpp>
#include <opengv/optimization_tools/solver_tools/SolverToolsNoncentralRelativePose.hpp>

#include "opencv2/core.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/xfeatures2d.hpp"

#include <opengv/amm.hpp>

#include <opengv/statistic/iterations_info.hpp>
#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"


using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;
void readme();
/*
 * @function main
 * @brief Main function
 */
int main( int argc, char** argv )
{
  if( argc != 3 )
  { readme(); return -1; }
  Mat img_1 = imread( argv[1], IMREAD_GRAYSCALE );
  Mat img_2 = imread( argv[2], IMREAD_GRAYSCALE );
  if( !img_1.data || !img_2.data )
  { std::cout<< " --(!) Error reading images " << std::endl; return -1; }
  //-- Step 1: Detect the keypoints using SURF Detector, compute the descriptors
  int minHessian = 400;
  Ptr<SURF> detector = SURF::create();
  detector->setHessianThreshold(minHessian);
  std::vector<KeyPoint> keypoints_1, keypoints_2;
  Mat descriptors_1, descriptors_2;
  detector->detectAndCompute( img_1, Mat(), keypoints_1, descriptors_1 );
  detector->detectAndCompute( img_2, Mat(), keypoints_2, descriptors_2 );
  //-- Step 2: Matching descriptor vectors using FLANN matcher
  FlannBasedMatcher matcher;
  std::vector< DMatch > matches;
  matcher.match( descriptors_1, descriptors_2, matches );
  double max_dist = 0; double min_dist = 100;
  //-- Quick calculation of max and min distances between keypoints
  for( int i = 0; i < descriptors_1.rows; i++ )
  { double dist = matches[i].distance;
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
  { if( matches[i].distance <= max(2*min_dist, 0.02) )
    { good_matches.push_back( matches[i]); }
  }
  //-- Draw only "good" matches
  /*Mat img_matches;
  drawMatches( img_1, keypoints_1, img_2, keypoints_2,
               good_matches, img_matches, Scalar::all(-1), Scalar::all(-1),
               vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
  //-- Show detected matches
  imshow( "Good Matches", img_matches );*/
  std::vector<KeyPoint> k1;
  std::vector<KeyPoint> k2;
  for( int i = 0; i < (int)good_matches.size(); i++ )
  {
    //printf( "-- Good Match [%d] Keypoint 1: %d  -- Keypoint 2: %d  \n", i, good_matches[i].queryIdx, good_matches[i].trainIdx );
    std::cout << ( "-- Good Match [%d] Keypoint 1: %d  -- Keypoint 2: %d  \n", i, good_matches[i].queryIdx, good_matches[i].trainIdx ) << std::endl;
    k1.push_back(keypoints_1[good_matches[i].queryIdx]);
    k2.push_back(keypoints_2[good_matches[i].trainIdx]);    
  }
  for( int i = 0; i < (int)good_matches.size(); i++ )
    {
      std::cout << "Point1: " << k1[i].pt << std::endl;
      std::cout << "Point2: " << k2[i].pt << std::endl;
    }

  std::cout << "Cal matrices: " << std::endl;
  Eigen::Matrix3d K_left = Eigen::Matrix3d::Identity(3,3);
  K_left(0,0) = 9.597910e+02; K_left(0,1) =  0.000000e+00; K_left(0,2) = 6.960217e+02;
  K_left(1,0) = 0.000000e+00; K_left(1,1) =  9.569251e+02; K_left(1,2) = 2.241806e+02;
  K_left(2,0) = 0.000000e+00; K_left(2,1) =  0.000000e+00; K_left(2,2) = 1.000000e+00;

  Eigen::Matrix3d K_right = Eigen::Matrix3d::Identity(3,3);
  K_right(0,0) = 9.037596e+02; K_right(0,1) =  0.000000e+00; K_right(0,2) = 6.957519e+02;
  K_right(1,0) = 0.000000e+00; K_right(1,1) =  9.019653e+02; K_right(1,2) = 2.242509e+02;
  K_right(2,0) = 0.000000e+00; K_right(2,1) =  0.000000e+00; K_right(2,2) = 1.000000e+00;

  //To use the methods it will be necessary to build the central adapter
  opengv::bearingVectors_t v1;
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
  opengv::relative_pose::CentralRelativeAdapter adapter(v1, v2);
  //timer
  struct timeval tic;
  struct timeval toc;
  opengv::rotation_t rot = Eigen::MatrixXd::Identity(3,3);
  opengv::translation_t trans = Eigen::MatrixXd::Zero(3,1);
  trans(0,0) = 0.0;trans(1,0) = 0.5;trans(2,0) = 0.0;
  ObjectiveFunctionInfo * obj_func = new SquaredFunctionNoIterationsInfo(adapter);

  SolverTools * solver_func = new SolverToolsNoncentralRelativePose();
  double tol = 1e-6;
  double step = 0.01;
  amm optimizer_func;
  std::vector<iterations_info> stat_iter;
  opengv::transformation_t result = optimizer_func.amm_solver(tol, rot, trans, obj_func, solver_func, step, stat_iter);
  std::cout << "Rotation: " << std::endl << result.block<3,3>(0,0) << std::endl;
  std::cout << "Translation: " << std::endl << result.block<3,1>(0,3) << std::endl;
  delete solver_func;
  delete obj_func;

  std::cout << "running fivept_nister" << std::endl;
  opengv::essentials_t fivept_nister_essentials;
  gettimeofday( &tic, 0 );
  fivept_nister_essentials = opengv::relative_pose::fivept_nister(adapter);
  gettimeofday( &toc, 0 );
  double fivept_nister_time = TIMETODOUBLE(timeval_minus(toc,tic));

   std::cout << "running fivept_stewenius" << std::endl;
   opengv::complexEssentials_t fivept_stewenius_essentials;
   gettimeofday( &tic, 0 );
   fivept_stewenius_essentials = opengv::relative_pose::fivept_stewenius(adapter);
   gettimeofday( &toc, 0 );
   double fivept_stewenius_time = TIMETODOUBLE(timeval_minus(toc,tic));
   std::cout << "results from stewenius' five-point algorithm:" << std::endl;
   for(int i = 0; i < (int) fivept_stewenius_essentials.size(); ++i){
     std::cout << fivept_stewenius_essentials.at(i) << std::endl << std::endl;
   }
   std::cout << "results from nisters' five-point algorithm:" << std::endl;
   for( size_t i = 0; i < fivept_nister_essentials.size(); i++ )
     std::cout << fivept_nister_essentials.at(i) << std::endl << std::endl;
   //waitKey(0);
   opengv::translation_t twopt_translation;
   gettimeofday( &tic, 0 );
   twopt_translation = opengv::relative_pose::twopt(adapter,true);
   gettimeofday( &toc, 0 );
   double twopt_time = TIMETODOUBLE(timeval_minus(toc,tic));
   std::cout << "results from two-points algorithm:" << std::endl;
   std::cout << twopt_translation << std::endl << std::endl;

   std::cout << "running eightpt" << std::endl;
   opengv::essential_t eightpt_essential;
   gettimeofday( &tic, 0 );
   eightpt_essential = opengv::relative_pose::eightpt(adapter);
   gettimeofday( &toc, 0 );
   double eightpt_time = TIMETODOUBLE(timeval_minus(toc,tic));
   std::cout << "results from eight-point algorithm:" << std::endl;
   std::cout << eightpt_essential << std::endl << std::endl;

   /*NON LINEAR*/
   std::cout << "setting perturbed rotation and ";
   std::cout << "running eigensolver" << std::endl;
   opengv::translation_t t_perturbed;
   opengv::rotation_t R_perturbed;
   Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity(3,3);
   Eigen::MatrixXd position = Eigen::MatrixXd::Zero(3,1);
   position(0,0) = 0.5;position(1,0) = 0.5;position(2,0) = 0.5;
   opengv::getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.01);
   opengv::rotation_t eigensolver_rotation;
   gettimeofday( &tic, 0 );
   adapter.setR12(R_perturbed);
   eigensolver_rotation = opengv::relative_pose::eigensolver(adapter);
   gettimeofday( &toc, 0 );
   double eigensolver_time = TIMETODOUBLE(timeval_minus(toc,tic));

   std::cout << "setting perturbed pose and ";
   std::cout << "performing nonlinear optimization" << std::endl;
   opengv::getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1);
   opengv::transformation_t nonlinear_transformation;
   gettimeofday( &tic, 0 );
   adapter.sett12(t_perturbed);
   adapter.setR12(R_perturbed);
   nonlinear_transformation = opengv::relative_pose::optimize_nonlinear(adapter);
   gettimeofday( &toc, 0 );
   double nonlinear_time = TIMETODOUBLE(timeval_minus(toc,tic));

   for(int i=5; i < good_matches.size(); ++i){
     std::cout << "setting perturbed pose and ";
     std::cout << "performing nonlinear optimization with " << i << " indices" << std::endl;
     std::vector<int> indices10 = opengv::getNindices(i);
     opengv::getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1);
     adapter.sett12(t_perturbed);
     adapter.setR12(R_perturbed);
     opengv::transformation_t nonlinear_transformation_10 =
       opengv::relative_pose::optimize_nonlinear(adapter,indices10);
     std::cout << "results from nonlinear algorithm with only few correspondences:";
     std::cout << std::endl;
     std::cout << std::endl << nonlinear_transformation_10 << std::endl;
   }
   /************/
   std::cout << "results from eigensystem based rotation solver:" << std::endl;
   std::cout << eigensolver_rotation << std::endl << std::endl << std::endl;
   std::cout << "results from nonlinear algorithm:" << std::endl;
   std::cout << nonlinear_transformation << std::endl << std::endl;
   return 0;
}
/*
 * @function readme
 */
void readme()
{ std::cout << " Usage: ./SURF_FlannMatcher <img1> <img2>" << std::endl; }
