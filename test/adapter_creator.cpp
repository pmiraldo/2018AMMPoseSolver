#include "adapter_creator.hpp"
#include "experiment_helpers.hpp"
#include <Eigen/SVD>
#include <iostream>

void create_random_pose(opengv::rotation_t & rotation, opengv::translation_t & position, const double & noise, const double & outlierFraction, const size_t & numberPoints, const int & numberCameras,opengv::bearingVectors_t & bearingVectors, std::vector<int> & camCorrespondences, opengv::points_t & points, opengv::translations_t & camOffsets, opengv::rotations_t & camRotations){
  
  //create a random camera-system
  opengv::generateRandomCameraSystem( numberCameras, camOffsets, camRotations );
  
  //derive correspondences based on random point-cloud
  Eigen::MatrixXd gt(3,numberPoints);
  opengv::generateRandom2D3DCorrespondences(position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction, bearingVectors, points, camCorrespondences, gt );
  
  //print the experiment characteristics
  opengv::printExperimentCharacteristics(
				 position, rotation, noise, outlierFraction );
}

void choose_best_transformation(const opengv::transformations_t &solutions, const opengv::transformation_t & ref, opengv::transformation_t & transformation_method){
  double error = 1e7;
  int index = 0;
  for(int i = 0; i < solutions.size(); ++i){
    if(error > (ref - solutions[i] ).norm()){
      error = (ref - solutions[i] ).norm();
      index = i;
    }
  }
  transformation_method = solutions[index];
 }

void choose_best_essential_matrix(const opengv::essentials_t & solutions, const opengv::essential_t & ref, opengv::essential_t & essential_matrix){
  if(solutions.size() == 1){
    essential_matrix = solutions[0];
  }
  double error = 1e7;
  int index = 0;
  for(int i = 0; i < solutions.size(); ++i){
    if(error > (ref - solutions[i]).norm()){
      error = (ref - solutions[i]).norm();
      index = i;
    }
  }
  essential_matrix = solutions[index];
}


opengv::essential_t calculate_essential_matrix(const opengv::translation_t & position, const opengv::rotation_t & rotation)
{
  //E transforms vectors from vp 2 to 1: x_1^T * E * x_2 = 0
  //and E = (t)_skew*R
  Eigen::Matrix3d t_skew = Eigen::Matrix3d::Zero();
  t_skew(0,1) = -position[2];
  t_skew(0,2) = position[1];
  t_skew(1,0) = position[2];
  t_skew(1,2) = -position[0];
  t_skew(2,0) = -position[1];
  t_skew(2,1) = position[0];
  opengv::essential_t E = t_skew * rotation;
  std::cout << "the random essential matrix is:" << std::endl;
  std::cout << E << std::endl;
  return E;
}


opengv::transformation_t calculate_transformation_from_essential_matrix(const opengv::essential_t & essential_matrix, const opengv::rotation_t &ref_rotation, const opengv::translation_t & ref_translation){
  //Obtain matrices U and V
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(essential_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix3d U = svd.matrixU();
  Eigen::Matrix3d V = svd.matrixV();
  Eigen::Matrix3d sing_values = Eigen::Matrix3d::Zero(3,3);
  sing_values(0,0) = svd.singularValues()(0,0);
  sing_values(1,1) = svd.singularValues()(1,0);

  double detU = U.determinant();
  if (detU < 0){
    U = -U;
  }
  
  double detV = V.determinant();
  if (detV < 0){
    V = -V;
  }
  
  Eigen::Matrix3d Rz = Eigen::Matrix3d::Zero(3,3);
  Rz(0,1) = -1;Rz(1,0) = 1; Rz(2,2) = 1;

  
  opengv::rotation_t r = U * Rz.transpose() * V.transpose();
  Eigen::MatrixXd skew_t = U * Rz * sing_values * U.transpose();
  Eigen::MatrixXd t = Eigen::MatrixXd::Zero(3,1);
  t(0,0) = skew_t(2,1);t(1,0) = skew_t(0,2); t(2,0) = skew_t(1,0);
  
  if ( (r - ref_rotation).norm() > 1e-6){
    r = U * Rz * V.transpose();
  }
  if( (t - ref_translation).norm() > 1e-6){
    skew_t = U * Rz.transpose() * sing_values * U.transpose();
    t(0,0) = skew_t(2,1);t(1,0) = skew_t(0,2); t(2,0) = skew_t(1,0);
  }
  
  opengv::transformation_t result;
  result.block<3,3>(0,0) = r;
  result.block<3,1>(0,3) = t;
  return result;
}
