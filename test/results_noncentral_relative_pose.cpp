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

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"

#include <opengv/optimization_tools/solver_tools/SolverToolsNoncentralRelativePose.hpp>
#include <opengv/optimization_tools/objective_function_tools/SquaredFunctionInfo.hpp>
#include <opengv/optimization_tools/objective_function_tools/SquaredFunctionNoIterationsInfo.hpp>

#include <opengv/statistic/AggregateStatisticalInfo.hpp>
#include <opengv/statistic/StatisticalInfoContainer.hpp>
#include <opengv/statistic/iterations_info.hpp>
#include <opengv/amm.hpp>



using namespace std;
using namespace Eigen;
using namespace opengv;

int main( int argc, char** argv )
{
  // initialize random seed
  initializeRandomSeed();
  int n_experiments = 30;//10;
  double percentage = 0.2;
  int noise_levels = 7;//4;

  //set experiment parameters
  double noise = 0.0;
  double outlierFraction = 0.0;
  size_t numberPoints = 50;
  int numberCameras = 8;
  std::ofstream error_file("relative_pose_error.csv");
  std::ofstream iterations_file("relative_pose_iterations.csv");
  std::vector<StatisticalInfoContainer> statistical_error_info_ge;
  std::vector<StatisticalInfoContainer> statistical_error_info_17pt;
  std::vector<StatisticalInfoContainer> statistical_error_info_amm;
  std::vector<StatisticalInfoContainer> statistical_error_info_nonlin;
  for(int index = 0; index < noise_levels; index++){
    //set noise
    double noise = 0.0 + 1 * index;
    std::cout << std::endl << std::endl << "***************************" << std::endl;
    std::cout << "Noise: " << noise << std::endl;
    int index_stat = 0;

    //This part will be used to estimate the average so we can eliminate outliers
    int total_realizations = 0;
    std::vector<double> average_estimator_rotation;
    std::vector<double> average_estimator_translation;
    int estimator_stat_counter = 0;
    while(estimator_stat_counter < percentage * n_experiments){
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
      SolverTools * solver_container = NULL;
      solver_container = new SolverToolsNoncentralRelativePose();
      ObjectiveFunctionInfo * info_container_noiterations = NULL;
      info_container_noiterations = new SquaredFunctionNoIterationsInfo(adapter);

      transformation_t seventeenpt_transformation_all;
      seventeenpt_transformation_all = relative_pose::seventeenpt(adapter);
      double tol = 1e-6;
      double step = 0.01;
      amm solver_object;
      std::vector<iterations_info> stat_iterations_list;
      transformation_t amm_solution_noiterations = solver_object.amm_solver( tol, seventeenpt_transformation_all.block<3,3>(0,0),
									     seventeenpt_transformation_all.block<3,1>(0,3),
									          info_container_noiterations,
									     solver_container, step, stat_iterations_list);

      double rotation_error_amm = (amm_solution_noiterations.block<3,3>(0,0) - rotation).norm();
      double translation_error_amm = (amm_solution_noiterations.block<3,1>(0,3) - position).norm();
      average_estimator_rotation.push_back(rotation_error_amm);
      average_estimator_translation.push_back(translation_error_amm);
      
      estimator_stat_counter++;
    }
    std::sort(average_estimator_rotation.rbegin(), average_estimator_rotation.rend());
    std::sort(average_estimator_translation.rbegin(), average_estimator_translation.rend());

    //Obtain the median
    double median_rotation = 0;
    double median_translation = 0;
    if (average_estimator_rotation.size() % 2 == 0){
      int index_1 = average_estimator_rotation.size() / 2;
      int index_2 = index_1 - 1;
      median_rotation = 0.5 * (average_estimator_rotation[index_1] + average_estimator_rotation[index_2]);
      median_translation = 0.5 * (average_estimator_translation[index_1] + average_estimator_translation[index_2]);
    }
    else{
      int index = (average_estimator_rotation.size() - 1) / 2;
      median_rotation = average_estimator_rotation[index];
      median_translation = average_estimator_translation[index];
    }
    for(int i = 0; i < average_estimator_rotation.size(); ++i){
      average_estimator_rotation[i] = average_estimator_rotation[i] - median_rotation;
      average_estimator_translation[i] = average_estimator_translation[i] - median_translation;
    }
    median_rotation = 0;
    median_translation = 0;
    if (average_estimator_rotation.size() % 2 == 0){
      int index_1 = average_estimator_rotation.size() / 2;
      int index_2 = index_1 - 1;
      median_rotation = 0.5 * (average_estimator_rotation[index_1] + average_estimator_rotation[index_2]);
      median_translation = 0.5 * (average_estimator_translation[index_1] + average_estimator_translation[index_2]);
    }
    else{
      int index = (average_estimator_rotation.size() - 1) / 2;
       median_rotation = average_estimator_rotation[index];
       median_translation = average_estimator_translation[index];
    }
    double MAD_rotation = 1.4826 * median_rotation;
    double MAD_translation = 1.4826 * median_translation;
    
    
      
    //End of the part that eliminates outliers
    while(index_stat < n_experiments ){
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
      for( size_t i = 0; i < iterations; i++ ){
	sixpt_rotations = relative_pose::sixpt(adapter,indices6);
      }
      gettimeofday( &toc, 0 );
      double sixpt_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
      std::cout << "running ge with 8 correspondences" << std::endl;
      std::vector<int> indices8 = getNindices(8);
      transformation_t ge_transformation;
      gettimeofday( &tic, 0 );
      for( size_t i = 0; i < iterations; i++ ){
	ge_transformation = relative_pose::ge(adapter,indices8);
      }
      gettimeofday( &toc, 0 );
      double ge_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
      std::cout << "running seventeenpt algorithm with 17 correspondences";
      std::cout << std::endl;
      std::vector<int> indices17 = getNindices(17);
      transformation_t seventeenpt_transformation;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++){
	seventeenpt_transformation = relative_pose::seventeenpt(adapter,indices17);
      }
      gettimeofday( &toc, 0 );
      double seventeenpt_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

      std::cout << "running seventeenpt algorithm with all correspondences";
      std::cout << std::endl;
      transformation_t seventeenpt_transformation_all;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++){
	seventeenpt_transformation_all = relative_pose::seventeenpt(adapter);
      }
      gettimeofday( &toc, 0 );
      double seventeenpt_time_all = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
      std::cout << "setting perturbed pose and ";
      std::cout << "performing nonlinear optimization" << std::endl;
      translation_t t_perturbed; rotation_t R_perturbed;
      getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1);
      transformation_t nonlinear_transformation;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++){
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
      transformation_t nonlinear_transformation_10 = relative_pose::optimize_nonlinear(adapter,indices10);
      
      //print results
      /*std::cout << "results from 6pt algorithm:" << std::endl;
      for( size_t i = 0; i < sixpt_rotations.size(); i++ ){
	std::cout << sixpt_rotations[i] << std::endl << std::endl;
      }
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
      std::cout << "AMM" << std::endl;*/

      //Create point with info
      /*ObjectiveFunctionInfo * info_container = NULL;
      info_container = new SquaredFunctionInfo(adapter);

      //Create solver pointer
      SolverTools * solver_container = NULL;
      solver_container = new SolverToolsNoncentralRelativePose();
      amm solver_object;
      gettimeofday(&tic,0);
      double step = 0.01;
      transformation_t amm_solution = solver_object.amm_solver( tol, initial_rotation, initial_translation, info_container, solver_container, step);
      gettimeofday(&toc, 0);
      double time_amm_solution = TIMETODOUBLE(timeval_minus(toc,tic)); //  / iterations;*/

      //NEW APPROACE SAME THING
      std::vector<iterations_info> stat_iterations_list;
      //Create point with info
      SolverTools * solver_container = NULL;
      solver_container = new SolverToolsNoncentralRelativePose();
      ObjectiveFunctionInfo * info_container_noiterations = NULL;
      info_container_noiterations = new SquaredFunctionNoIterationsInfo(adapter);
      
      double tol = 1e-6;
      double step = 0.01;
      amm solver_object;
      gettimeofday(&tic,0);
      transformation_t amm_solution_noiterations = solver_object.amm_solver( tol, seventeenpt_transformation_all.block<3,3>(0,0),
									     seventeenpt_transformation_all.block<3,1>(0,3),
									          info_container_noiterations,
									     solver_container, step, stat_iterations_list);
      gettimeofday(&toc, 0);
      //delete info_container;
      delete solver_container;
      delete info_container_noiterations;

      
      double time_amm_solution_noiterations = TIMETODOUBLE(timeval_minus(toc,tic)); //  / iterations;
      double rotation_error_amm = (amm_solution_noiterations.block<3,3>(0,0) - rotation).norm();
      double translation_error_amm = (amm_solution_noiterations.block<3,1>(0,3) - position).norm();
      StatisticalInfoContainer trial_statistical_info_amm(noise, "amm", rotation_error_amm, translation_error_amm, time_amm_solution_noiterations, stat_iterations_list);
      
      stat_iterations_list.clear();
      double rotation_error_17pt = (seventeenpt_transformation_all.block<3,3>(0,0) - rotation).norm();
      double translation_error_17pt = (seventeenpt_transformation_all.block<3,1>(0,3) - position).norm();
      StatisticalInfoContainer trial_statistical_info_17pt(noise, "17 pt", rotation_error_17pt, translation_error_17pt, seventeenpt_time_all, stat_iterations_list);

      double rotation_error_ge = (ge_transformation.block<3,3>(0,0) - rotation).norm();
      double translation_error_ge = (ge_transformation.block<3,1>(0,3) - position).norm();
      StatisticalInfoContainer trial_statistical_info_ge(noise, "ge", rotation_error_ge, translation_error_ge, ge_time, stat_iterations_list);

      double rotation_error_nonlin = (nonlinear_transformation.block<3,3>(0,0) - rotation).norm();
      double translation_error_nonlin = (nonlinear_transformation.block<3,1>(0,3) - position).norm();
      StatisticalInfoContainer trial_statistical_info_nonlin(noise, "nonlin", rotation_error_nonlin, translation_error_nonlin, nonlinear_time, stat_iterations_list);

      statistical_error_info_ge.push_back(trial_statistical_info_ge);
      statistical_error_info_17pt.push_back(trial_statistical_info_17pt);
      statistical_error_info_amm.push_back(trial_statistical_info_amm);
      statistical_error_info_nonlin.push_back(trial_statistical_info_nonlin);
      index_stat++;
    }
    total_realizations++;
  }
  for(int i = 0; i < statistical_error_info_ge.size(); ++i){
    statistical_error_info_ge[i].printInfo(error_file, iterations_file, false);
  }
  for(int i = 0; i < statistical_error_info_17pt.size(); ++i){
    statistical_error_info_17pt[i].printInfo(error_file, iterations_file, false);
  }
  for(int i = 0; i < statistical_error_info_nonlin.size(); ++i){
    statistical_error_info_nonlin[i].printInfo(error_file, iterations_file, false);
  }
  for(int i = 0; i < statistical_error_info_amm.size(); ++i){
    statistical_error_info_amm[i].printInfo(error_file, iterations_file, true);
  }
  iterations_file.close();
  error_file.close();
}
