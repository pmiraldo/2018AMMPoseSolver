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
  int n_experiments = 150;//10;
   int noise_levels = 10;//4;

  //set experiment parameters
  double noise = 0.0;
  double outlierFraction = 0.0;
  size_t numberPoints = 100;
  int numberCameras = 4;
  std::ofstream error_file("relative_pose_error.csv");
  std::ofstream iterations_file("relative_pose_iterations.csv");
  std::vector<std::vector<StatisticalInfoContainer> > statistical_error_methods(4);
  for(int index = 0; index < noise_levels; index++){
    //set noise
    double noise = 0.0 + 1 * index;
    std::cout << std::endl << std::endl << "***************************" << std::endl;
    std::cout << "Noise: " << noise << std::endl;
    int index_stat = 0;

    //This part will be used to estimate the average so we can eliminate outliers
    int total_realizations = 0;
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

       transformation_t gt_transformation;
      gt_transformation.block<3,3>(0,0) = rotation;
      gt_transformation.block<3,1>(0,3) = position;

      //timer
      struct timeval tic;
      struct timeval toc;
      size_t iterations = 10;

      std::cout << "running ge with all correspondences" << std::endl;
      transformation_t ge_transformation;
      gettimeofday( &tic, 0 );
      for( size_t i = 0; i < iterations; i++ )
	ge_transformation = relative_pose::ge(adapter);
      gettimeofday( &toc, 0 );
      double ge_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
     
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
      transformation_t amm_solution_noiterations;
      gettimeofday(&tic,0);
      for(int i = 0; i < iterations; ++i){
	amm_solution_noiterations = solver_object.amm_solver( tol, seventeenpt_transformation_all.block<3,3>(0,0),
							      seventeenpt_transformation_all.block<3,1>(0,3),
							      info_container_noiterations,
							      solver_container, step, stat_iterations_list);
      }
      gettimeofday(&toc, 0);
      //delete info_container;
      delete solver_container;
      delete info_container_noiterations;
      double time_amm_solution_noiterations = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      
      StatisticalInfoContainer trial_statistical_info_amm(noise, "amm", amm_solution_noiterations, gt_transformation, time_amm_solution_noiterations, stat_iterations_list);
      stat_iterations_list.clear();
      StatisticalInfoContainer trial_statistical_info_17pt(noise, "17 pt", seventeenpt_transformation_all, gt_transformation, seventeenpt_time_all, stat_iterations_list);
      StatisticalInfoContainer trial_statistical_info_ge(noise, "ge", ge_transformation, gt_transformation, ge_time, stat_iterations_list);
      StatisticalInfoContainer trial_statistical_info_nonlin(noise, "nonlin", nonlinear_transformation, gt_transformation, nonlinear_time, stat_iterations_list);

      statistical_error_methods[0].push_back(trial_statistical_info_amm);
      statistical_error_methods[1].push_back(trial_statistical_info_17pt);
      statistical_error_methods[2].push_back(trial_statistical_info_ge);
      statistical_error_methods[3].push_back(trial_statistical_info_nonlin);
      index_stat++;
    }
    total_realizations++;
  }
  /*for(int i = 0; i < statistical_error_info_ge.size(); ++i){
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
    }*/
  iterations_file.close();
  error_file.close();
}
