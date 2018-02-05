#include "adapter_creator.hpp"
#include "experiment_helpers.hpp"

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

void choose_best_rotation(const opengv::transformations_t &solutions, const opengv::rotation_t & rotation, opengv::transformation_t & transformation_method){
  double error = 1e7;
  int index = 0;
  for(int i = 0; i < solutions.size(); ++i){
    if(error > (rotation - solutions[i].block<3,3>(0,0) ).norm()){
      error = (rotation - solutions[i].block<3,3>(0,0) ).norm();
      index = i;
    }
  }
  transformation_method = solutions[index];
}
