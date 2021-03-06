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

#include <cmath>
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
#include <random>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"

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
  int n_experiments = 2;//500;
  int noise_levels = 3;//4;
  std::vector<std::string> methods_list;
  methods_list.push_back("gp3p");
  methods_list.push_back("gpnp");
  methods_list.push_back("upnp");
  methods_list.push_back("non linear");
  methods_list.push_back("amm gpnp");
  methods_list.push_back("amm upnp");
  methods_list.push_back("amm inf norm");
  int position_table_stat = 0;
  StatisticalInfoContainer statistical_info(n_experiments, noise_levels, 2, methods_list);
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

      std::cout << "THE OFFSET IS: " << position_table_stat << std::endl;
      std::vector<iterations_info> iterations_list;
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
      choose_best_transformation(gp3p_transformations, gt_transformation, gp3p_transformation);
      statistical_info.append_info(0, position_table_stat, noise, gp3p_transformation, gt_transformation, gp3p_time, iterations_list);
      
      std::cout << "running gpnp over all correspondences" << std::endl;
      transformation_t gpnp_transformation;
      gettimeofday( &tic, 0 );
      for(size_t i = 0; i < iterations; i++){
	gpnp_transformation = absolute_pose::gpnp(adapter);
      }
      gettimeofday( &toc, 0 );
      double gpnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;    
      statistical_info.append_info(1, position_table_stat, noise, gpnp_transformation, gt_transformation, gpnp_time, iterations_list);
      
      std::cout << "running upnp over all correspondences" << std::endl;
      transformations_t upnp_transformations;
      gettimeofday( &tic, 0 );
      for( size_t i = 0; i < iterations; i++ ){
	upnp_transformations = absolute_pose::upnp(adapter);
      }
      gettimeofday( &toc, 0);
      double upnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      transformation_t upnp_transformation;
      choose_best_transformation(upnp_transformations, gt_transformation, upnp_transformation);
      statistical_info.append_info(2, position_table_stat, noise, upnp_transformation, gt_transformation, upnp_time, iterations_list);
      
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
      statistical_info.append_info(3, position_table_stat, noise, nonlinear_transformation, gt_transformation, nonlinear_time, iterations_list);
      
      
      
      //Create solver pointer and solver tools
      SolverTools * solver_container = solver_container = new SolverToolsNoncentralRelativePose();
      amm solver_object;
      std::default_random_engine generator;
      std::normal_distribution<double> distribution(0.0, 9 * (index + 1));

      double angle = (M_PI / 20 ) * distribution(generator);
      std::cout << "The angle is: " << angle << std::endl;
      rotation_t error_rot = Eigen::Matrix3d::Identity(3,3);
      error_rot(0,0) = std::cos(angle);error_rot(0,1) = -std::sin(angle);
      error_rot(1,0) = std::sin(angle);error_rot(1,1) = std::cos(angle);
      std::cout << "Rotation error: " << std::endl << error_rot << std::endl;
      rotation_t rotation_init = rotation * error_rot;
      translation_t translation_init = position;
      //AMM GlobalPnPFunctionInfo
      double tol = 1e-6;
      double step = 0.008;
      transformation_t global_pnp;
      ObjectiveFunctionInfo * info_container_amm_gpnp = new GlobalPnPFunctionInfo(adapter);
      gettimeofday(&tic,0);
      for(int i = 0; i < iterations; ++i){
	global_pnp = solver_object.amm_solver( tol, rotation_init.inverse(), -rotation_init.inverse() * translation_init, info_container_amm_gpnp, solver_container, step, iterations_list);
      }
      gettimeofday(&toc,0);
      delete info_container_amm_gpnp;
      double time_pnp_global = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      statistical_info.append_info(4, position_table_stat, noise, global_pnp, gt_transformation, time_pnp_global, iterations_list);

      //AMM OptimalPnPFunctionInfo
      step = 0.0051;
      transformation_t optimal_upnp;
      
      ObjectiveFunctionInfo * info_container_optimal_upnp = new OptimalUPnPFunctionInfo(adapter);
      gettimeofday(&tic, 0);
      for(int i = 0; i < iterations; ++i){
	optimal_upnp = solver_object.amm_solver( tol, rotation_init.inverse(), -rotation_init.inverse() * translation_init, info_container_optimal_upnp, solver_container, step, iterations_list);
      }
      gettimeofday(&toc, 0);
      delete info_container_optimal_upnp;
      double time_optimal_upnp = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      statistical_info.append_info(5, position_table_stat, noise, optimal_upnp, gt_transformation, time_optimal_upnp, iterations_list);
      
      //AMM GlobalPnPFunction Infinite norm
      step = 0.45;
      tol = 1e-6;
      transformation_t infinite_norm_gpnp;
      ObjectiveFunctionInfo * info_container_inf_norm = new GlobalPnPInfiniteNormFunctionInfo(adapter, rotation_init.inverse(), -rotation_init.inverse() * translation_init);
      gettimeofday(&tic, 0);
      for(int i = 0; i < iterations; ++i){
	infinite_norm_gpnp = solver_object.amm_solver( tol, rotation_init.inverse(), -rotation_init.inverse() * translation_init, info_container_inf_norm, solver_container, step, iterations_list);
      }
      gettimeofday(&toc, 0);
      delete info_container_inf_norm;
      double time_infinite_norm_gpnp = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
      statistical_info.append_info(6, position_table_stat, noise, infinite_norm_gpnp, gt_transformation, time_infinite_norm_gpnp, iterations_list);
      
      //The solver is no longer needed. So it can be erased
      delete solver_container;
      iterations_list.clear();
      std::cout << "Errors: "        << std::endl;
      std::cout << "gp3p (rot) "         << (gp3p_transformation.block<3,3>(0,0) - gt_transformation.block<3,3>(0,0) ).norm() << std::endl;
      std::cout << "gpnp (rot) "         << (gpnp_transformation.block<3,3>(0,0) - gt_transformation.block<3,3>(0,0) ).norm() << std::endl;
      std::cout << "upnp (rot) "         << (upnp_transformation.block<3,3>(0,0) - gt_transformation.block<3,3>(0,0) ).norm() << std::endl;
      std::cout << "nonlin (rot) "       << (nonlinear_transformation.block<3,3>(0,0) - gt_transformation.block<3,3>(0,0) ).norm() << std::endl;
      std::cout << "amm gpnp (rot) "     << (global_pnp.block<3,3>(0,0) - gt_transformation.block<3,3>(0,0) ).norm() << std::endl;
      std::cout << "amm upnp (rot) "     << (optimal_upnp.block<3,3>(0,0) - gt_transformation.block<3,3>(0,0) ).norm() << std::endl;
      std::cout << "amm inf norm (rot) " << (infinite_norm_gpnp.block<3,3>(0,0) - gt_transformation.block<3,3>(0,0) ).norm() << std::endl;

      std::cout << "Errors: "        << std::endl;
      std::cout << "gp3p (trans) "         << (gp3p_transformation.block<3,1>(0,3)               - gt_transformation.block<3,1>(0,3) ).norm() << std::endl;
      std::cout << "gpnp (trans) "         << (gpnp_transformation.block<3,1>(0,3)               - gt_transformation.block<3,1>(0,3) ).norm() << std::endl;
      std::cout << "upnp (trans) "         << (upnp_transformation.block<3,1>(0,3)               - gt_transformation.block<3,1>(0,3) ).norm() << std::endl;
      std::cout << "nonlin (trans) "       << (nonlinear_transformation.block<3,1>(0,3)          - gt_transformation.block<3,1>(0,3) ).norm() << std::endl;
      std::cout << "amm gpnp (trans) "     << (global_pnp.block<3,1>(0,3)                        - gt_transformation.block<3,1>(0,3) ).norm() << std::endl;
      std::cout << "amm upnp (trans) "     << (optimal_upnp.block<3,1>(0,3)                      - gt_transformation.block<3,1>(0,3) ).norm() << std::endl;
      std::cout << "amm inf norm (trans) " << (infinite_norm_gpnp.block<3,1>(0,3) - gt_transformation.block<3,1>(0,3) ).norm() << std::endl;

      std::cout << "Time: "       << std::endl;
      std::cout << "gp3p"         <<  gp3p_time         << std::endl;
      std::cout << "gpnp"         <<  gpnp_time         << std::endl;
      std::cout << "upnp"         <<  upnp_time         << std::endl;
      std::cout << "nonlin"       <<  nonlinear_time    << std::endl;
      std::cout << "amm gpnp"     <<  time_pnp_global   << std::endl;
      std::cout << "amm upnp"     <<  time_optimal_upnp << std::endl;
      std::cout << "amm inf norm" <<  time_infinite_norm_gpnp << std::endl;
      
      index_stat++;
      position_table_stat++;
    }
  }
  std::cout << "Write in file" << std::endl;
  std::string filename = "info.txt";
  statistical_info.print_info(filename);
}

