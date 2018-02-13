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
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <sstream>
#include <fstream>

#include "adapter_creator.hpp"
#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"

#include <opengv/optimization_tools/solver_tools/SolverToolsNoncentralRelativePose.hpp>
#include <opengv/optimization_tools/objective_function_tools/SquaredFunctionCentralInfo.hpp>

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
  int n_experiments = 500;
  int noise_levels = 10;
  
  //set experiment parameters
  double outlierFraction = 0.0;
  size_t numberPoints = 100;
  std::ofstream error_file("relative_pose_error.csv");
  std::ofstream iterations_file("relative_pose_iterations.csv");
  std::vector<std::vector<StatisticalInfoContainer> > statistical_error_methods(6);
  for(int index = 0; index < noise_levels; index++){
    //set noise
    double noise = 0.0 + 1 * index;
    std::cout << std::endl << std::endl << "***************************" << std::endl;
    std::cout << "Noise: " << noise << std::endl;
    int index_stat = 0;
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
       generateCentralCameraSystem( camOffsets, camRotations );
    
       //derive correspondences based on random point-cloud
       bearingVectors_t bearingVectors1;
       bearingVectors_t bearingVectors2;
       std::vector<int> camCorrespondences1; //unused in the central case
       std::vector<int> camCorrespondences2; //unused in the central case
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
       
       //compute and print the essential-matrix
       //printEssentialMatrix( position, rotation );
       essential_t essential_matrix_ref = calculate_essential_matrix(position, rotation);
       //create a central relative adapter
       relative_pose::CentralRelativeAdapter adapter(
						  bearingVectors1,
						  bearingVectors2,
						  rotation);

       transformation_t gt_transformation = calculate_transformation_from_essential_matrix(essential_matrix_ref, rotation, position);
       
       //timer
       struct timeval tic;
       struct timeval toc;
       size_t iterations = 10;

       
       std::vector<iterations_info> stat_iterations_list;
       //running experiments
       std::cout << "running twopt" << std::endl;
       translation_t twopt_translation;
       gettimeofday( &tic, 0 );
       for(size_t i = 0; i < iterations; i++)
	 twopt_translation = relative_pose::twopt(adapter,true);
       gettimeofday( &toc, 0 );
       double twopt_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
       transformation_t twopt_transformation = Eigen::MatrixXd::Zero(3,4);
       twopt_transformation.block<3,1>(0,3) = twopt_translation;
       StatisticalInfoContainer trial_statistical_info_twopt(noise, "2pt", twopt_transformation, gt_transformation, twopt_time, stat_iterations_list);

       
       std::cout << "running sevenpt" << std::endl;
       essentials_t sevenpt_essentials;
       essential_t sevenpt_essential;
       gettimeofday( &tic, 0 );
       for(size_t i = 0; i < iterations; i++)
	 sevenpt_essentials = relative_pose::sevenpt(adapter);
       gettimeofday( &toc, 0 );
       double sevenpt_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
       sevenpt_essential = sevenpt_essentials[0];
       choose_best_essential_matrix(sevenpt_essentials, essential_matrix_ref, sevenpt_essential);
       transformation_t sevenpt_transformation = calculate_transformation_from_essential_matrix(sevenpt_essential, rotation, position);
       StatisticalInfoContainer trial_statistical_info_sevenpt(noise, "7pt", sevenpt_transformation, gt_transformation, sevenpt_time, stat_iterations_list);
       
       std::cout << "running eightpt" << std::endl;
       essential_t eightpt_essential;
       gettimeofday( &tic, 0 );
       for(size_t i = 0; i < iterations; i++)
	 eightpt_essential = relative_pose::eightpt(adapter);
       gettimeofday( &toc, 0 );
       double eightpt_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
       transformation_t eightpt_transformation = calculate_transformation_from_essential_matrix(eightpt_essential, rotation, position);
       StatisticalInfoContainer trial_statistical_info_eightpt(noise, "8pt", eightpt_transformation, gt_transformation, eightpt_time, stat_iterations_list);
       
       std::cout << "setting perturbed rotation and ";
       std::cout << "running eigensolver" << std::endl;
       translation_t t_perturbed; rotation_t R_perturbed;
       getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.01);
       rotation_t eigensolver_rotation;
       gettimeofday( &tic, 0 );
       for(size_t i = 0; i < iterations; i++)
	 {
	   adapter.setR12(R_perturbed);
	   eigensolver_rotation = relative_pose::eigensolver(adapter);
	 }
       gettimeofday( &toc, 0 );
       double eigensolver_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
       transformation_t eigen_transformation = Eigen::MatrixXd::Zero(3,4);
       eigen_transformation.block<3,3>(0,0) = eigensolver_rotation;
       StatisticalInfoContainer trial_statistical_info_eigen(noise, "eigen", eigen_transformation, gt_transformation, eigensolver_time, stat_iterations_list);
       
       std::cout << "setting perturbed pose and ";
       std::cout << "performing nonlinear optimization" << std::endl;
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
       StatisticalInfoContainer trial_statistical_info_nonlin(noise, "nonlin", nonlinear_transformation, gt_transformation, nonlinear_time, stat_iterations_list);
       
       //NEW APPROACE SAME THING
       
       //Create point with info
       SolverTools * solver_container = NULL;
       solver_container = new SolverToolsNoncentralRelativePose();
       ObjectiveFunctionInfo * info_container = NULL;
       info_container = new SquaredFunctionCentralInfo(adapter);
       
       double tol = 1e-12;
       double step = 0.05;
       amm solver_object;
       transformation_t amm_solution;
       rotation_t rotation_init = R_perturbed;
       translation_t translation_init = t_perturbed;
       Eigen::MatrixXd trans = twopt_translation;
       gettimeofday(&tic,0);
       for(int i = 0; i < iterations; ++i){
	 amm_solution = solver_object.amm_solver( tol, R_perturbed, t_perturbed, info_container,  solver_container, step, stat_iterations_list);
       }
       gettimeofday(&toc, 0);
       delete info_container;
       delete solver_container;
       double time_amm_solution = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

       StatisticalInfoContainer trial_statistical_info_amm(noise, "amm", amm_solution, gt_transformation, time_amm_solution, stat_iterations_list);

       //print results
       /*std::cout << "results from seven-point algorithm:" << std::endl;
       for (int i = 0; i < sevenpt_essentials.size(); ++i){
	 std::cout << "solution: " << i << std::endl;
	 std::cout << std::endl << sevenpt_essentials[i] << std::endl;
       }
       
       
       std::cout << "results from two-points algorithm:" << std::endl;
       std::cout << twopt_translation << std::endl << std::endl;
       std::cout << "Solution chosen for sevenpt: " << std::endl;
       std::cout << "Essential matrix 7 pt" << std::endl;
       for(int i = 0; i < sevenpt_essentials.size(); ++i){
	 std::cout << sevenpt_essentials[i] << std::endl;
       }
       std::cout << "7 pt transformation" << std::endl;
       std::cout << sevenpt_transformation << std::endl;
       std::cout << "gt transformation" << std::endl;
       std::cout << gt_transformation << std::endl;
       std::cout << "results from eight-point algorithm:" << std::endl;
       std::cout << eightpt_transformation << std::endl << std::endl;
       std::cout << "results from eigensystem based rotation solver:" << std::endl;
       std::cout << eigensolver_rotation << std::endl << std::endl << std::endl;
       std::cout << "results from nonlinear algorithm:" << std::endl;
       std::cout << nonlinear_transformation << std::endl;
       std::cout << "results from amm: " << std::endl;
       std::cout << amm_solution << std::endl << std::endl*/;

       /*std::cout << "error rot   amm " << (gt_transformation.block<3,3>(0,0) - amm_solution.block<3,3>(0,0)).norm() << std::endl;
       std::cout << "error rot   7pt " << (gt_transformation.block<3,3>(0,0) - sevenpt_transformation.block<3,3>(0,0)).norm() << std::endl;
       std::cout << "error rot   8pt " << (gt_transformation.block<3,3>(0,0) - eightpt_transformation.block<3,3>(0,0)).norm() << std::endl;
       std::cout << "error rot eigen " << (gt_transformation.block<3,3>(0,0) - eigen_transformation.block<3,3>(0,0)).norm() << std::endl;
       std::cout << "error rot nolin " << (gt_transformation.block<3,3>(0,0) - nonlinear_transformation.block<3,3>(0,0)).norm() << std::endl;

       std::cout << "error trabs   2pt " << (gt_transformation.block<3,1>(0,3) - twopt_transformation.block<3,1>(0,3)).norm() << std::endl;
       std::cout << "error trabs   amm " << (gt_transformation.block<3,1>(0,3) - amm_solution.block<3,1>(0,3)).norm() << std::endl;
       std::cout << "error trans   7pt " << (gt_transformation.block<3,1>(0,3) - sevenpt_transformation.block<3,1>(0,3)).norm() << std::endl;
       std::cout << "error trans   8pt " << (gt_transformation.block<3,1>(0,3) - eightpt_transformation.block<3,1>(0,3)).norm() << std::endl;
       std::cout << "error trans eigen " << (gt_transformation.block<3,1>(0,3) - eigen_transformation.block<3,1>(0,3)).norm() << std::endl;
       std::cout << "error trans nolin " << (gt_transformation.block<3,1>(0,3) - nonlinear_transformation.block<3,1>(0,3)).norm() << std::endl;*/
       
       /*std::cout << "timings from seven-point algorithm: ";
       std::cout << sevenpt_time << std::endl;
       std::cout << "timings from eight-point algorithm: ";
       std::cout << eightpt_time << std::endl;
       std::cout << "timings from eigensystem based rotation solver: ";
       std::cout << eigensolver_time << std::endl;
       std::cout << "timings from nonlinear algorithm: ";
       std::cout << nonlinear_time << std::endl;
       std::cout << "timing amm" << std::endl;
       std::cout << time_amm_solution << std::endl;*/
       
       statistical_error_methods[0].push_back(trial_statistical_info_twopt);
       statistical_error_methods[1].push_back(trial_statistical_info_sevenpt);
       statistical_error_methods[2].push_back(trial_statistical_info_eightpt);
       statistical_error_methods[3].push_back(trial_statistical_info_eigen);
       statistical_error_methods[4].push_back(trial_statistical_info_nonlin);
       statistical_error_methods[5].push_back(trial_statistical_info_amm);
       
       index_stat++;
     }
  }
  bool flag = false;
  for(int i = 0; i < statistical_error_methods.size(); ++i){
    if (i == statistical_error_methods.size() - i){
      flag = true;
    }
    for(int j = 0; j < statistical_error_methods[i].size(); ++j){     
      statistical_error_methods[i][j].printInfo(error_file, iterations_file, flag);
    }
  }
  iterations_file.close();
  error_file.close();
}
