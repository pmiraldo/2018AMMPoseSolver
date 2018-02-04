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
#include <opengv/amm.hpp>


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
  int n_experiments = 20;//10;
  int noise_levels = 3;//4;
  std::ofstream error_file("absolute_pose_error.csv");
  std::ofstream iterations_file("absolute_pose_iterations.csv");
  std::vector<StatisticalInfoContainer> statistical_error_info_gp3p;
  std::vector<StatisticalInfoContainer> statistical_error_info_gpnp;
  std::vector<StatisticalInfoContainer> statistical_error_info_upnp;
  std::vector<StatisticalInfoContainer> statistical_error_info_nonlin;
  std::vector<StatisticalInfoContainer> statistical_error_info_amm_global_pnp;
  std::vector<StatisticalInfoContainer> statistical_error_info_amm_optimal_upnp;
  std::vector<StatisticalInfoContainer> statistical_error_info_amm_infinite_norm;
 
  for(int index = 0; index < noise_levels; index++){

    double noise = 0.0 + 1 * index;
    int index_stat = 0;
    int total_realizations = 0;
    while(index_stat < n_experiments){
      //create a random viewpoint pose
      translation_t position = generateRandomTranslation(2.0);
      rotation_t rotation = generateRandomRotation(0.5);
  
      //create a random camera-system
      translations_t camOffsets;
      rotations_t camRotations;
      generateRandomCameraSystem( numberCameras, camOffsets, camRotations );
  
      //derive correspondences based on random point-cloud
      bearingVectors_t bearingVectors;
      points_t points;
      std::vector<int> camCorrespondences;
      Eigen::MatrixXd gt(3,numberPoints);
      generateRandom2D3DCorrespondences(
					position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction,
					bearingVectors, points, camCorrespondences, gt );

      //print the experiment characteristics
      printExperimentCharacteristics(
				     position, rotation, noise, outlierFraction );

      //create a non-central absolute adapter
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
      size_t iterations = 50;

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

      std::cout << "running gpnp over all correspondences" << std::endl;
      transformation_t gpnp_transformation;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++){
	gpnp_transformation = absolute_pose::gpnp(adapter);
      }
      gettimeofday( &toc, 0 );
      double gpnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

      std::cout << "running gpnp over 6 correspondences" << std::endl;
      std::vector<int> indices6 = getNindices(6);
      transformation_t gpnp_transformation_6 =
	absolute_pose::gpnp( adapter, indices6 );

      std::cout << "running upnp over all correspondences" << std::endl;
      transformations_t upnp_transformations;
      gettimeofday( &tic, 0 );
      for( size_t i = 0; i < iterations; i++ ){
	upnp_transformations = absolute_pose::upnp(adapter);
      }
      gettimeofday( &toc, 0);
      double upnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
      std::cout << "running upnp over 3 correspondences" << std::endl;
      std::vector<int> indices3 = getNindices(3);
      transformations_t upnp_transformations_3 =
	absolute_pose::upnp( adapter, indices3 );

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

      std::cout << "setting perturbed pose and ";
      std::cout << "performing nonlinear optimization with 10 correspondences";
      std::cout << std::endl;
      std::vector<int> indices10 = getNindices(10);
      //add a small perturbation to the rotation
      getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1 );
      adapter.sett(t_perturbed);
      adapter.setR(R_perturbed);
      transformation_t nonlinear_transformation_10 =
	absolute_pose::optimize_nonlinear(adapter,indices10);

      //print the results
      std::cout << "results from gp3p algorithm:" << std::endl;
      for(size_t i = 0; i < gp3p_transformations.size(); i++){
	std::cout << gp3p_transformations[i] << std::endl << std::endl;
      }
      std::cout << "results from gpnp algorithm:" << std::endl;
      std::cout << gpnp_transformation << std::endl << std::endl;
      std::cout << "results from gpnp algorithm with only 6 correspondences:";
      std::cout << std::endl;
      std::cout << gpnp_transformation_6 << std::endl << std::endl;
      std::cout << "results from upnp algorithm:" << std::endl;
      for(size_t i = 0; i < upnp_transformations.size(); i++){
	std::cout << upnp_transformations[i] << std::endl << std::endl;
      }
      std::cout << "results from upnp algorithm with only 3 correspondences:";
      std::cout << std::endl;
      for(size_t i = 0; i < upnp_transformations_3.size(); i++){
	std::cout << upnp_transformations_3[i] << std::endl << std::endl;
      }
      std::cout << "results from nonlinear algorithm:" << std::endl;
      std::cout << nonlinear_transformation << std::endl << std::endl;
      std::cout << "results from nonlinear algorithm with only 10 correspondences:";
      std::cout << std::endl;
      std::cout << nonlinear_transformation_10 << std::endl << std::endl;

      /*std::cout << "timings from gp3p algorithm: ";
      std::cout << gp3p_time << std::endl;
      std::cout << "timings from gpnp algorithm: ";
      std::cout << gpnp_time << std::endl;
      std::cout << "timings from upnp algorithm: ";
      std::cout << upnp_time << std::endl;
      std::cout << "timings from nonlinear algorithm: ";
      std::cout << nonlinear_time << std::endl;*/

      std::cout << "Verification Info" << std::endl;
      //Choose the right initial state for 
      rotation_t rotation_gp3p;
      translation_t translation_gp3p;
      double error = 1e7;
      int index = 0;
      //Choose the right solution from the gp3p method
      for(int i = 0; i < gp3p_transformations.size(); ++i){
	if(error < (rotation - gp3p_transformations[i].block<3,3>(0,0) ).norm()){
	  error = (rotation - gp3p_transformations[i].block<3,3>(0,0) ).norm();
	  index = i;
	}
      }
      rotation_gp3p    = gp3p_transformations[index].block<3,3>(0,0);
      translation_gp3p = gp3p_transformations[index].block<3,1>(0,0);
      
      
      std::vector<iterations_info> iterations_list;
      //Create solver pointer and solver tools
      SolverTools * solver_container = solver_container = new SolverToolsNoncentralRelativePose();
      amm solver_object;

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

      double rotation_error_amm_global_pnp = (global_pnp.block<3,3>(0,0) - rotation).norm();
      double translation_error_amm_global_pnp = (global_pnp.block<3,1>(0,3) - position).norm();
      StatisticalInfoContainer trial_statistical_info_amm_global_pnp(noise, "amm gpnp", rotation_error_amm_global_pnp, translation_error_amm_global_pnp, time_pnp_global, iterations_list);


      
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
      double rotation_error_amm_optimal_upnp = (optimal_upnp.block<3,3>(0,0) - rotation).norm();
      double translation_error_amm_optimal_upnp = (optimal_upnp.block<3,1>(0,3) - position).norm();
      StatisticalInfoContainer trial_statistical_info_amm_optimal_upnp(noise, "amm upnp", rotation_error_amm_global_pnp, translation_error_amm_optimal_upnp, time_optimal_upnp, iterations_list);

      //AMM GlobalPnPFunction Infinite norm
      step = 0.45;
      tol = 1e-15;
      transformation_t infinite_norm_gpnp;
      gettimeofday(&tic, 0);
     
      for(int i = 0; i < iterations; ++i){
	ObjectiveFunctionInfo * info_container = new GlobalPnPInfiniteNormFunctionInfo(adapter, rotation_gp3p.inverse(), -rotation_gp3p.inverse() * translation_gp3p);
	infinite_norm_gpnp = solver_object.amm_solver( tol, rotation_gp3p.inverse(), -rotation_gp3p.inverse() * translation_gp3p, info_container, solver_container, step, iterations_list);
	delete info_container;
      }
      gettimeofday(&toc, 0);
      double time_infinite_norm_gpnp = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

      double rotation_error_amm_infinite_norm_gpnp = (infinite_norm_gpnp.block<3,3>(0,0) - rotation).norm();
      double translation_error_amm_infinite_norm_gpnp = (infinite_norm_gpnp.block<3,1>(0,3) - position).norm();
      StatisticalInfoContainer trial_statistical_info_amm_infinite_norm_gpnp(noise, "amm infinite norm", rotation_error_amm_infinite_norm_gpnp, translation_error_amm_infinite_norm_gpnp, time_infinite_norm_gpnp, iterations_list);
      //The solver is no longer needed. So it can be erased
      delete solver_container;


      //Choose the right solution for the upnp 
      rotation_t rotation_upnp;
      translation_t translation_upnp;
      error = 1e7;
      index = 0;
      //Choose the right solution from the upnp method
      for(int i = 0; i < upnp_transformations.size(); ++i){
	if(error < (rotation - upnp_transformations[i].block<3,3>(0,0) ).norm()){
	  error = (rotation - upnp_transformations[i].block<3,3>(0,0) ).norm();
	  index = i;
	}
      }
      rotation_upnp    = upnp_transformations[index].block<3,3>(0,0);
      translation_upnp = upnp_transformations[index].block<3,1>(0,0);

      iterations_list.clear();
      index_stat++;
      statistical_error_info_amm_global_pnp.push_back(trial_statistical_info_amm_global_pnp);
      statistical_error_info_amm_optimal_upnp.push_back(trial_statistical_info_amm_optimal_upnp);
      statistical_error_info_amm_infinite_norm.push_back(trial_statistical_info_amm_infinite_norm_gpnp);

      double rotation_error_gp3p = (rotation_gp3p - rotation).norm();
      double translation_error_gp3p = (translation_gp3p - position).norm();
      double rotation_error_upnp = (rotation_upnp - rotation).norm();
      double translation_error_upnp = (translation_upnp - position).norm();
      double rotation_error_gpnp = (gpnp_transformation.block<3,3>(0,0) - rotation).norm();
      double translation_error_gpnp = (gpnp_transformation.block<3,1>(0,3) - position).norm();
      double rotation_error_nonlin = (nonlinear_transformation.block<3,3>(0,0) - rotation).norm();
      double translation_error_nonlin = (nonlinear_transformation.block<3,1>(0,3) - position).norm();

      StatisticalInfoContainer trial_statistical_info_gp3p(noise, "gp3p", rotation_error_gp3p, translation_error_gp3p, gp3p_time, iterations_list);
      StatisticalInfoContainer trial_statistical_info_upnp(noise, "upnp", rotation_error_upnp, translation_error_upnp, upnp_time, iterations_list);
      StatisticalInfoContainer trial_statistical_info_gpnp(noise, "gpnp", rotation_error_gpnp, translation_error_gpnp, gpnp_time, iterations_list);
      StatisticalInfoContainer trial_statistical_info_nonlin(noise, "nonlin", rotation_error_nonlin, translation_error_nonlin, nonlinear_time, iterations_list);
      statistical_error_info_gp3p.push_back(trial_statistical_info_gp3p);
      statistical_error_info_gpnp.push_back(trial_statistical_info_upnp);
      statistical_error_info_upnp.push_back(trial_statistical_info_gpnp);
      statistical_error_info_nonlin.push_back(trial_statistical_info_nonlin);
    }
    total_realizations++;
  }
  for(int i = 0; i < statistical_error_info_amm_global_pnp.size(); ++i){
    statistical_error_info_amm_global_pnp[i].printInfo(error_file, iterations_file, true);
  }
  for(int i = 0; i < statistical_error_info_amm_optimal_upnp.size(); ++i){
    statistical_error_info_amm_optimal_upnp[i].printInfo(error_file, iterations_file, true);
  }
  for(int i = 0; i < statistical_error_info_amm_infinite_norm.size(); ++i){
    statistical_error_info_amm_infinite_norm[i].printInfo(error_file, iterations_file, true);
  }
  for(int i = 0; i < statistical_error_info_gp3p.size(); ++i){
    statistical_error_info_gp3p[i].printInfo(error_file, iterations_file, false);
  }
  for(int i = 0; i < statistical_error_info_gpnp.size(); ++i){
    statistical_error_info_gpnp[i].printInfo(error_file, iterations_file, false);
  }
  for(int i = 0; i < statistical_error_info_upnp.size(); ++i){
    statistical_error_info_upnp[i].printInfo(error_file, iterations_file, false);
  }
  for(int i = 0; i < statistical_error_info_nonlin.size(); ++i){
    statistical_error_info_nonlin[i].printInfo(error_file, iterations_file, false);
  }
  iterations_file.close();
  error_file.close();
}

