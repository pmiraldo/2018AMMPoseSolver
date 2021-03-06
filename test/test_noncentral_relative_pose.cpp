/******************************************************************************
 * Author:   Laurent Kneip                                                    *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/NoncentralRelativeAdapter.hpp>
#include <sstream>
#include <fstream>
#include <math.h>
#include <random>

#include <opengv/amm.hpp>


#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"


using namespace std;
using namespace Eigen;
using namespace opengv;

int main( int argc, char** argv )
{
  // initialize random seed
  initializeRandomSeed();

  //set experiment parameters
  double noise = 0.0;
  double outlierFraction = 0.0;
  size_t numberPoints = 100;
  int numberCameras = 4;

  //generate a random pose for viewpoint 1
  translation_t position1 = Eigen::Vector3d::Zero();
  rotation_t rotation1 = Eigen::Matrix3d::Identity();

  //generate a random pose for viewpoint 2
  translation_t position2 = generateRandomTranslation(2.0);
  rotation_t rotation2 = generateRandomRotation(0.5);

  //create a fake central camera
  translations_t camOffsets;
  rotations_t camRotations;
  generateRandomCameraSystem( numberCameras, camOffsets, camRotations );

  //derive correspondences based on random point-cloud
  bearingVectors_t bearingVectors1;
  bearingVectors_t bearingVectors2;
  std::vector<int> camCorrespondences1;
  std::vector<int> camCorrespondences2;
  Eigen::MatrixXd gt(3,numberPoints);
  generateRandom2D2DCorrespondences(
      position1, rotation1, position2, rotation2,
      camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors1, bearingVectors2,
      camCorrespondences1, camCorrespondences2, gt );

  //Extract the relative pose
  translation_t position; rotation_t rotation;
  extractRelativePose(
      position1, position2, rotation1, rotation2, position, rotation, false );

  //print experiment characteristics
  printExperimentCharacteristics( position, rotation, noise, outlierFraction );

  //create non-central relative adapter
  relative_pose::NoncentralRelativeAdapter adapter(
      bearingVectors1,
      bearingVectors2,
      camCorrespondences1,
      camCorrespondences2,
      camOffsets,
      camRotations,
      position,
      rotation);

  //timer
  struct timeval tic;
  struct timeval toc;
  size_t iterations = 100;

  //running experiment
  std::cout << "running sixpt with 6 correspondences" << std::endl;
  std::vector<int> indices6 = getNindices(6);
  rotations_t sixpt_rotations;
  gettimeofday( &tic, 0 );
  for( size_t i = 0; i < iterations; i++ )
    sixpt_rotations = relative_pose::sixpt(adapter,indices6);
  gettimeofday( &toc, 0 );
  double sixpt_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
  std::cout << "running ge with 8 correspondences" << std::endl;
  std::vector<int> indices8 = getNindices(8);
  transformation_t ge_transformation;
  gettimeofday( &tic, 0 );
  for( size_t i = 0; i < iterations; i++ )
    ge_transformation = relative_pose::ge(adapter,indices8);
  gettimeofday( &toc, 0 );
  double ge_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
  std::cout << "running seventeenpt algorithm with 17 correspondences";
  std::cout << std::endl;
  std::vector<int> indices17 = getNindices(17);
  transformation_t seventeenpt_transformation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    seventeenpt_transformation = relative_pose::seventeenpt(adapter,indices17);
  gettimeofday( &toc, 0 );
  double seventeenpt_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running seventeenpt algorithm with all correspondences";
  std::cout << std::endl;
  transformation_t seventeenpt_transformation_all;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    seventeenpt_transformation_all = relative_pose::seventeenpt(adapter);
  gettimeofday( &toc, 0 );
  double seventeenpt_time_all =
      TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
  std::cout << "setting perturbed pose and ";
  std::cout << "performing nonlinear optimization" << std::endl;
  translation_t t_perturbed; rotation_t R_perturbed;
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1);
  transformation_t nonlinear_transformation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
  {
    adapter.sett12(t_perturbed);
    adapter.setR12(R_perturbed);
    nonlinear_transformation = relative_pose::optimize_nonlinear(adapter);
  }
  gettimeofday( &toc, 0 );
  double nonlinear_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "setting perturbed pose and ";
  std::cout << "performing nonlinear optimization with 10 correspondences";
  std::cout << std::endl;
  std::vector<int> indices10 = getNindices(10);
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1);
  adapter.sett12(t_perturbed);
  adapter.setR12(R_perturbed);
  transformation_t nonlinear_transformation_10 =
      relative_pose::optimize_nonlinear(adapter,indices10);

  //print results
  std::cout << "results from 6pt algorithm:" << std::endl;
  for( size_t i = 0; i < sixpt_rotations.size(); i++ )
    std::cout << sixpt_rotations[i] << std::endl << std::endl;
  std::cout << "result from ge using 8 points:" << std::endl;
  std::cout << ge_transformation << std::endl << std::endl;
  std::cout << "results from 17pt algorithm:" << std::endl;
  std::cout << seventeenpt_transformation << std::endl << std::endl;
  std::cout << "results from 17pt algorithm with all points:" << std::endl;
  std::cout << seventeenpt_transformation_all << std::endl << std::endl;
  std::cout << "results from nonlinear algorithm:" << std::endl;
  std::cout << nonlinear_transformation << std::endl << std::endl;
  std::cout << "results from nonlinear algorithm with only few correspondences:";
  std::cout << std::endl;
  std::cout << nonlinear_transformation_10 << std::endl << std::endl;
  
  std::cout << "timings from 6pt algorithm: ";
  std::cout << sixpt_time << std::endl;
  std::cout << "timings from ge: ";
  std::cout << ge_time << std::endl;
  std::cout << "timings from 17pt algorithm: ";
  std::cout << seventeenpt_time << std::endl;
  std::cout << "timings from 17pt algorithm with all the points: ";
  std::cout << seventeenpt_time_all << std::endl;
  std::cout << "timings from nonlinear algorithm: ";
  std::cout << nonlinear_time << std::endl;
  std::cout << "AMM" << std::endl;
  double tol = 1e-12;
  rotation_t initial_state = MatrixXd::Identity(3,3);

  /*transformation_t amm_solution;
 
  std::mt19937 rng;
  rng.seed(std::random_device()());
  std::uniform_int_distribution<std::mt19937::result_type> dist(0.01,1); // distribution in range [0, 1]
  double angle =  M_PI / 36;
  double wx = 0; double wy = 0; double wz = 1;
  Eigen::Matrix3d skew_matrix = Eigen::Matrix3d::Zero(3,3);
  skew_matrix(0,1) = -wz; skew_matrix(0,2) = wy;
  skew_matrix(1,0) =  wz; skew_matrix(1,2) = -wx;
  skew_matrix(2,0) = -wy; skew_matrix(2,1) =  wx;
  rotation_t error_rotation = Eigen::Matrix3d::Identity(3,3) +
                              (std::sin(angle) / angle)*skew_matrix +
                              ( ( 1 - std::cos(angle) * std::cos(angle) ) / (angle * angle) ) * skew_matrix * skew_matrix;

  translation_t error_translation = Eigen::Vector3d::Random(3,1);
  error_translation = error_translation / error_translation.norm();
  rotation_t  initial_rotation = rotation * error_rotation;
  translation_t initial_translation = position;// + error_translation;
  gettimeofday(&tic,0);
  amm solver;
  double step = 0.1;
  for(int i = 0; i < iterations; i++){
    amm_solution  = solver.amm_solver(tol, initial_rotation, initial_translation, step);//initial_state);
  };
  gettimeofday(&toc,01); 
  double time_amm_solution =
  TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;*/

  /*std::cout << "********************************************************" << std::endl;
  std::cout << "Vector error: " << std::endl << skew_matrix << std::endl;
  std::cout << "Angle:        " << std::endl << angle       << std::endl;
  std::cout << "Erro da rotação: "    << std::endl << error_rotation    << std::endl;
  std::cout << "Erro da translação: " << std::endl << error_translation << std::endl;
  std::cout << "Solution presented by algorithm amm: " << std::endl;
  std::cout << amm_solution << std::endl;
  std::cout << "Real rotation : " << std::endl << rotation            << std::endl;
  std::cout << "Corrupted rotation: " << std::endl << initial_rotation    << std::endl;
  std::cout << "Initial translation: " << std::endl << initial_translation << std::endl;
  std::cout << "Error AMM: " << std::endl;
  std::cout << (amm_solution.block<3,3>(0,0) - rotation).norm() << std::endl;
  std::cout << "Time: " << time_amm_solution << std::endl;*/

  std::cout << "Error 17: " << std::endl;
  std::cout << (seventeenpt_transformation.block<3,3>(0,0) - rotation).norm() << std::endl;
  std::cout << "Time: " << seventeenpt_time_all << std::endl;
  std::cout << "Error nonlin: " << std::endl;
  std::cout << (nonlinear_transformation.block<3,3>(0,0) - rotation).norm() << std::endl;
  std::cout << "Time: " << nonlinear_time << std::endl;
  std::cout << "GE: " << std::endl;
  std::cout << (ge_transformation.block<3,3>(0,0) - rotation).norm() << std::endl;
  std::cout << "Time: " << ge_time << std::endl;

  std::cout << "GE Error translation: " << std::endl;
  std::cout << (ge_transformation.block<3,1>(0,3) - position).norm() << std::endl;
  std::cout << "Error 17 translation: " << std::endl;
  std::cout << (seventeenpt_transformation.block<3,1>(0,3) - position).norm() << std::endl;
  std::cout << "Error nonlin translation: " << std::endl;
  std::cout << (nonlinear_transformation.block<3,1>(0,3) - position).norm() << std::endl;
}
