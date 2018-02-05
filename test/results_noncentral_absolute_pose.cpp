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
#include <opengv/amm.hpp>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/NoncentralAbsoluteAdapter.hpp>
#include <opengv/optimization_tools/objective_function_tools/GlobalPnPFunctionInfo.hpp>
#include <opengv/optimization_tools/objective_function_tools/OptimalUPnPFunctionInfo.hpp>
#include <opengv/optimization_tools/objective_function_tools/GlobalPnPInfiniteNormFunctionInfo.hpp>
#include <opengv/optimization_tools/solver_tools/SolverToolsNoncentralRelativePose.hpp>
#include <sstream>
#include <fstream>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"

#include <opengv/statistic/AggregateStatisticalInfo.hpp>
#include <opengv/statistic/StatisticalInfoContainer.hpp>
#include <opengv/statistic/iterations_info.hpp>

#include "adapter_creator.hpp"
using namespace std;
using namespace Eigen;
using namespace opengv;

int main( int argc, char** argv )
{
  //initialize random seed
  initializeRandomSeed();
  
  //set experiment parameters
  double noise = 0.0;
  double outlierFraction = 0.0;
  size_t numberPoints = 100;
  int numberCameras = 4;

  //Experience parameters
  int n_experiments = 150;
  int noise_levels = 10;//4;
  std::ofstream error_file("absolute_pose_error.csv");
  std::ofstream iterations_file("absolute_pose_iterations.csv");
  std::vector<std::vector<StatisticalInfoContainer> > statistical_error_methods(7);
  for(int index = 0; index < noise_levels; index++){

    double noise = 0.0 + 1 * index;
    int index_stat = 0;
    int total_realizations = 0;
    
    while(index_stat < n_experiments){
      
      //create a random viewpoint pose
      translation_t position = generateRandomTranslation(2.0);
      rotation_t rotation = generateRandomRotation(0.5);
      transformation_t gt_transformation;
      gt_transformation.block<3,3>(0,0) = rotation;
      gt_transformation.block<3,1>(0,3) = position;
      transformation_t gt_transformation_amm;
      gt_transformation_amm.block<3,3>(0,0) = rotation.inverse();
      gt_transformation_amm.block<3,1>(0,3) = -rotation.inverse() * position;
      
      bearingVectors_t bearingVectors;
      std::vector<int> camCorrespondences;
      points_t points;
      translations_t camOffsets;
      rotations_t camRotations;
      //create a non-central absolute adapter
      create_random_pose(rotation, position, noise, outlierFraction, numberPoints, numberCameras, bearingVectors, camCorrespondences, points, camOffsets, camRotations);
      absolute_pose::NoncentralAbsoluteAdapter adapter(
						       bearingVectors,
						       camCorrespondences,
						       points,
						       camOffsets,
						       camRotations,
						       rotation);
      
      //timer
      struct timeval tic;
      struct timeval toc;
      size_t iterations = 10;

      //run the experiments
      std::cout << "running Kneip's GP3P (using first three correspondences/";
      std::cout << std::endl;
      transformations_t gp3p_transformations;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++){
	gp3p_transformations = absolute_pose::gp3p(adapter);
      }
      gettimeofday( &toc, 0 );
      double gp3p_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      transformation_t gp3p_transformation;
      choose_best_rotation(gp3p_transformations, rotation, gp3p_transformation);
      

      std::cout << "running gpnp over all correspondences" << std::endl;
      transformation_t gpnp_transformation;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++){
	gpnp_transformation = absolute_pose::gpnp(adapter);
      }
      gettimeofday( &toc, 0 );
      double gpnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

      std::cout << "running upnp over all correspondences" << std::endl;
      transformations_t upnp_transformations;
      gettimeofday( &tic, 0 );
      for( size_t i = 0; i < iterations; i++ ){
	upnp_transformations = absolute_pose::upnp(adapter);
      }
      gettimeofday( &toc, 0);
      double upnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      transformation_t upnp_transformation;
      choose_best_rotation(upnp_transformations, rotation, upnp_transformation);
      
      std::cout << "setting perturbed pose and ";
      std::cout << "performing nonlinear optimization" << std::endl;
      //add a small perturbation to the rotation
      translation_t t_perturbed; rotation_t R_perturbed;
      getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1 );
      transformation_t nonlinear_transformation;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++)
	{
	  adapter.sett(t_perturbed);
	  adapter.setR(R_perturbed);
	  nonlinear_transformation = absolute_pose::optimize_nonlinear(adapter);
	}
      gettimeofday( &toc, 0 );
      double nonlinear_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

      
      
      std::vector<iterations_info> iterations_list;
      //Create solver pointer and solver tools
      SolverTools * solver_container = solver_container = new SolverToolsNoncentralRelativePose();
      amm solver_object;
      rotation_t rotation_gp3p = gp3p_transformation.block<3,3>(0,0);
      translation_t translation_gp3p = gp3p_transformation.block<3,1>(0,3);
      //AMM GlobalPnPFunctionInfo
      double tol = 1e-6;
      double step = 0.008;
      transformation_t global_pnp;
      gettimeofday(&tic,0);
      for(int i = 0; i < iterations; ++i){
	ObjectiveFunctionInfo * info_container = new GlobalPnPFunctionInfo(adapter);
	global_pnp = solver_object.amm_solver( tol, rotation_gp3p.inverse(), -rotation_gp3p.inverse() * translation_gp3p, info_container, solver_container, step, iterations_list);
	delete info_container;
      }
      gettimeofday(&toc,0);
      double time_pnp_global = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

      //AMM OptimalPnPFunctionInfo
      step = 0.0051;
      transformation_t optimal_upnp;
      gettimeofday(&tic, 0);
      for(int i = 0; i < iterations; ++i){
	ObjectiveFunctionInfo * info_container = new OptimalUPnPFunctionInfo(adapter);
	optimal_upnp = solver_object.amm_solver( tol, rotation_gp3p.inverse(), -rotation_gp3p.inverse() * translation_gp3p, info_container, solver_container, step, iterations_list);
	delete info_container;
      }
      gettimeofday(&toc, 0);
      double time_optimal_upnp = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
     

      //AMM GlobalPnPFunction Infinite norm
      step = 0.45;
      tol = 1e-6;
      transformation_t infinite_norm_gpnp;
      gettimeofday(&tic, 0);
     
      for(int i = 0; i < iterations; ++i){
	ObjectiveFunctionInfo * info_container = new GlobalPnPInfiniteNormFunctionInfo(adapter, rotation_gp3p.inverse(), -rotation_gp3p.inverse() * translation_gp3p);
	infinite_norm_gpnp = solver_object.amm_solver( tol, rotation_gp3p.inverse(), -rotation_gp3p.inverse() * translation_gp3p, info_container, solver_container, step, iterations_list);
	delete info_container;
      }
      gettimeofday(&toc, 0);
      double time_infinite_norm_gpnp = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

      
      //The solver is no longer needed. So it can be erased
      delete solver_container;
      iterations_list.clear();


      index_stat++;
      
      StatisticalInfoContainer trial_statistical_info_gp3p(noise, "gp3p", gp3p_transformation, gt_transformation, gp3p_time, iterations_list);
      StatisticalInfoContainer trial_statistical_info_upnp(noise, "upnp", upnp_transformation, gt_transformation, upnp_time, iterations_list);
      StatisticalInfoContainer trial_statistical_info_gpnp(noise, "gpnp", gpnp_transformation, gt_transformation, gpnp_time, iterations_list);
      StatisticalInfoContainer trial_statistical_info_nonlin(noise, "nonlin", nonlinear_transformation, gt_transformation, nonlinear_time, iterations_list);
      StatisticalInfoContainer trial_statistical_info_amm_global_pnp(noise, "amm gpnp", global_pnp, gt_transformation_amm, time_pnp_global, iterations_list);
      StatisticalInfoContainer trial_statistical_info_amm_optimal_upnp(noise, "amm upnp", optimal_upnp, gt_transformation_amm, time_optimal_upnp, iterations_list);
      StatisticalInfoContainer trial_statistical_info_amm_infinite_norm_gpnp(noise, "amm infinite norm", infinite_norm_gpnp, gt_transformation_amm, time_infinite_norm_gpnp, iterations_list);
     
     
      statistical_error_methods[0].push_back(trial_statistical_info_gp3p);
      statistical_error_methods[1].push_back(trial_statistical_info_upnp);
      statistical_error_methods[2].push_back(trial_statistical_info_gpnp);
      statistical_error_methods[3].push_back(trial_statistical_info_nonlin);
      statistical_error_methods[4].push_back(trial_statistical_info_amm_global_pnp);
      statistical_error_methods[5].push_back(trial_statistical_info_amm_optimal_upnp);
      statistical_error_methods[6].push_back(trial_statistical_info_amm_infinite_norm_gpnp);
      
    }
  }
  bool is_amm = false;
  for(int i = 0; i < statistical_error_methods.size(); ++i){
    if(i > 3){
      is_amm = true;
    }
    for(int j = 0; j < statistical_error_methods[i].size(); ++j){
    statistical_error_methods[i][j].printInfo(error_file, iterations_file, is_amm);
    }
  }
  iterations_file.close();
  error_file.close();
}

