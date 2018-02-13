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
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/math/cayley.hpp>
#include <sstream>
#include <fstream>

#include <opengv/amm.hpp>
#include <opengv/optimization_tools/objective_function_tools/GlobalPnPFunctionInfo.hpp>
#include <opengv/optimization_tools/objective_function_tools/OptimalUPnPFunctionInfo.hpp>
#include <opengv/optimization_tools/objective_function_tools/GlobalPnPInfiniteNormFunctionInfo.hpp>
#include <opengv/optimization_tools/solver_tools/SolverToolsNoncentralRelativePose.hpp>

#include "adapter_creator.hpp"
#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"

#include <opengv/statistic/StatisticalInfoContainer.hpp>
#include <opengv/statistic/iterations_info.hpp>


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

  
  //Experience parameters
  int n_experiments = 500;
  int noise_levels = 10;
  std::ofstream error_file("absolute_pose_central_error.csv");
  std::ofstream iterations_file("absolute_pose_central_iterations.csv");
  std::vector<std::vector<StatisticalInfoContainer> > statistical_error_methods(8);
  for(int index = 0; index < noise_levels; index++){
    double noise = 0.0 + 1 * index;
    int index_stat = 0;
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
  
      //create a fake central camera
      translations_t camOffsets;
      rotations_t camRotations;
      generateCentralCameraSystem( camOffsets, camRotations );
  
      //derive correspondences based on random point-cloud
      bearingVectors_t bearingVectors;
      points_t points;
      std::vector<int> camCorrespondences; //unused in the central case!
      Eigen::MatrixXd gt(3,numberPoints);
      generateRandom2D3DCorrespondences(
					position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction,
					bearingVectors, points, camCorrespondences, gt );

      //print the experiment characteristics
      printExperimentCharacteristics(
				     position, rotation, noise, outlierFraction );

      //create a central absolute adapter
      absolute_pose::CentralAbsoluteAdapter adapter(
						    bearingVectors,
						    points,
						    rotation );

      //timer
      struct timeval tic;
      struct timeval toc;
      size_t iterations = 10;

      std::vector<iterations_info> iterations_list;
      //run the experiments
      std::cout << "running Kneip's P3P (first three correspondences)" << std::endl;
      transformations_t p3p_kneip_transformations;
      transformation_t p3p_kneip_transformation;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++)
	p3p_kneip_transformations = absolute_pose::p3p_kneip(adapter);
      gettimeofday( &toc, 0 );
      double p3p_kneip_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      choose_best_transformation(p3p_kneip_transformations, gt_transformation, p3p_kneip_transformation);
      StatisticalInfoContainer trial_statistical_info_p3p_kneip(noise, "p3p kneip", p3p_kneip_transformation, gt_transformation, p3p_kneip_time, iterations_list);

      
      std::cout << "running Gao's P3P (first three correspondences)" << std::endl;
      transformations_t p3p_gao_transformations;
      transformation_t p3p_gao_transformation;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++)
	p3p_gao_transformations = absolute_pose::p3p_gao(adapter);
      gettimeofday( &toc, 0 );
      double p3p_gao_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      choose_best_transformation(p3p_gao_transformations, gt_transformation, p3p_gao_transformation);
      StatisticalInfoContainer trial_statistical_info_p3p_gao(noise, "p3p gao", p3p_gao_transformation, gt_transformation, p3p_gao_time, iterations_list);
      

      std::cout << "running epnp (all correspondences)" << std::endl;
      transformation_t epnp_transformation;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++)
	epnp_transformation = absolute_pose::epnp(adapter);
      gettimeofday( &toc, 0 );
      double epnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      StatisticalInfoContainer trial_statistical_info_epnp(noise, "epnp", epnp_transformation, gt_transformation, epnp_time, iterations_list);
      
      std::cout << "running upnp with all correspondences" << std::endl;
      transformations_t upnp_transformations;
      transformation_t upnp_transformation;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++)
	upnp_transformations = absolute_pose::upnp(adapter);
      gettimeofday( &toc, 0 );
      double upnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      choose_best_transformation(upnp_transformations, gt_transformation, upnp_transformation);
      StatisticalInfoContainer trial_statistical_info_upnp(noise, "upnp", upnp_transformation, gt_transformation, upnp_time, iterations_list);

      std::cout << "setting perturbed pose";
      std::cout << "and performing nonlinear optimization" << std::endl;
      //add a small perturbation to the pose
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
      StatisticalInfoContainer trial_statistical_info_nonlin(noise, "nonlin", nonlinear_transformation, gt_transformation, nonlinear_time, iterations_list);

      //Create solver pointer and solver tools
      SolverTools * solver_container = solver_container = new SolverToolsNoncentralRelativePose();
      amm solver_object;
      std::default_random_engine generator;
      std::normal_distribution<double> distribution(0.0, 9 * (index + 1));

      double angle = (M_PI / 90 ) * distribution(generator);
      std::cout << "The angle is: " << angle << std::endl;
      rotation_t error_rot = Eigen::Matrix3d::Identity(3,3);
      error_rot(0,0) = std::cos(angle);error_rot(0,1) = -std::sin(angle);
      error_rot(1,0) = std::sin(angle);error_rot(1,1) = std::cos(angle);
      std::cout << "Rotation error: " << std::endl << error_rot << std::endl;
      rotation_t rotation_init = rotation * error_rot;
      translation_t translation_init = position;
      rotation_init = rotation_init.inverse();//rotation_init.inverse();
      translation_init = -rotation_init * translation_init;
      
      //AMM GlobalPnPFunctionInfo
      double tol = 1e-9;
      double step = 0.008;
      transformation_t global_pnp;
      ObjectiveFunctionInfo * info_container_amm_gpnp = new GlobalPnPFunctionInfo(adapter);
      gettimeofday(&tic,0);
      for(int i = 0; i < iterations; ++i){
	global_pnp = solver_object.amm_solver( tol, rotation_init, translation_init, info_container_amm_gpnp, solver_container, step, iterations_list);
      }
      gettimeofday(&toc,0);
      delete info_container_amm_gpnp;
      double time_pnp_global = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      StatisticalInfoContainer trial_statistical_info_amm_global_pnp(noise, "amm gpnp", global_pnp, gt_transformation_amm, time_pnp_global, iterations_list);
      
      
      //AMM OptimalPnPFunctionInfo
      step = 0.0051;
      transformation_t optimal_upnp;
      ObjectiveFunctionInfo * info_container_optimal_upnp = new OptimalUPnPFunctionInfo(adapter);
      gettimeofday(&tic, 0);
      for(int i = 0; i < iterations; ++i){
	optimal_upnp = solver_object.amm_solver( tol, rotation_init, translation_init, info_container_optimal_upnp, solver_container, step, iterations_list);
      }
      gettimeofday(&toc, 0);
      delete info_container_optimal_upnp;
      double time_optimal_upnp = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      StatisticalInfoContainer trial_statistical_info_amm_optimal_upnp(noise, "amm upnp", optimal_upnp, gt_transformation_amm, time_optimal_upnp, iterations_list);
      
      //AMM GlobalPnPFunction Infinite norm
      step = 0.45;
      tol = 1e-6;
      transformation_t infinite_norm_gpnp;
      ObjectiveFunctionInfo * info_container_inf_norm = new GlobalPnPInfiniteNormFunctionInfo(adapter, rotation_init, translation_init);
      gettimeofday(&tic, 0);
      for(int i = 0; i < iterations; ++i){
	infinite_norm_gpnp = solver_object.amm_solver( tol, rotation_init.inverse(), -rotation_init.inverse() * translation_init, info_container_inf_norm, solver_container, step, iterations_list);
      }
      gettimeofday(&toc, 0);
      delete info_container_inf_norm;
      double time_infinite_norm_gpnp = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      StatisticalInfoContainer trial_statistical_info_amm_infinite_norm_gpnp(noise, "amm infinite norm", infinite_norm_gpnp, gt_transformation_amm, time_infinite_norm_gpnp, iterations_list);
      
      //The solver is no longer needed. So it can be erased
      delete solver_container;
      iterations_list.clear();
      
      statistical_error_methods[0].push_back(trial_statistical_info_p3p_kneip);
      statistical_error_methods[1].push_back(trial_statistical_info_p3p_gao);
      statistical_error_methods[2].push_back(trial_statistical_info_epnp);
      statistical_error_methods[3].push_back(trial_statistical_info_upnp);
      statistical_error_methods[4].push_back(trial_statistical_info_nonlin);
      statistical_error_methods[5].push_back(trial_statistical_info_amm_global_pnp);
      statistical_error_methods[6].push_back(trial_statistical_info_amm_optimal_upnp);
      statistical_error_methods[7].push_back(trial_statistical_info_amm_infinite_norm_gpnp);
      index_stat++;
      //print the results
      /*std::cout << "results from Kneip's P3P algorithm:" << std::endl;
      std::cout << p3p_kneip_transformation << std::endl << std::endl;
      std::cout << "results from Gao's P3P algorithm:" << std::endl;
      std::cout << p3p_gao_transformation << std::endl << std::endl;
      std::cout << "results from epnp algorithm:" << std::endl;
      std::cout << epnp_transformation << std::endl << std::endl;
      std::cout << "results from upnp:" << std::endl;
      std::cout << upnp_transformation << std::endl << std::endl;
      std::cout << "results from nonlinear algorithm:" << std::endl;
      std::cout << nonlinear_transformation << std::endl << std::endl;
      std::cout << "the real transformation: " << std::endl;
      std::cout << gt_transformation << std::endl;
      std::cout << "the real transformation (amm)" << std::endl;
      std::cout << gt_transformation_amm << std::endl;
      std::cout << "results from gpnp amm" << std::endl;
      std::cout << global_pnp << std::endl << std::endl;
      std::cout << "results from upnp amm" << std::endl;
      std::cout << optimal_upnp << std::endl << std::endl;
      std::cout << "results from infinite norm gpnp" << std::endl;
      std::cout << infinite_norm_gpnp << std::endl << std::endl;
      
      std::cout << "timings from Kneip's P3P algorithm: ";
      std::cout << p3p_kneip_time << std::endl;
      std::cout << "timings from Gao's P3P algorithm: ";
      std::cout << p3p_gao_time << std::endl;
      std::cout << "timings from epnp algorithm: ";
      std::cout << epnp_time << std::endl;
      std::cout << "timings for the upnp algorithm: ";
      std::cout << upnp_time << std::endl;
      std::cout << "timings from nonlinear algorithm: ";
      std::cout << nonlinear_time << std::endl;
      std::cout << "timing for gpnp amm" << std::endl;
      std::cout << time_pnp_global << std::endl;
      std::cout << "timing for upnp amm" << std::endl;
      std::cout << time_optimal_upnp << std::endl;
      std::cout << "timing for infinite norm gpnp amm" << std::endl;
      std::cout << time_infinite_norm_gpnp << std::endl;*/
    }
  }
  bool is_amm = false;
  for(int i = 0; i < statistical_error_methods.size(); ++i){
    if(i > 4){
      is_amm = true;
    }
    for(int j = 0; j < statistical_error_methods[i].size(); ++j){
      statistical_error_methods[i][j].printInfo(error_file, iterations_file, is_amm);
    }
  }
  iterations_file.close();
  error_file.close();
} 
