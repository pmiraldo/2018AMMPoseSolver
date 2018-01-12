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
#include <opengv/optimization_tools/solver_tools/SolverToolsNoncentralRelativePose.hpp>
#include <sstream>
#include <fstream>

#include "container.h"
#include "statistic_info.h"
#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"


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
  int n_experiments = 1;//10;
  int noise_levels = 1;//4;

  std::vector<Container> information_gp3p;
  std::vector<Container> information_gpnp;
  std::vector<Container> information_upnp;
  std::vector<Container> information_nonlinear;
  std::vector<Container> information_amm;
  std::vector<Statistic_info> information_statistics;
  for(int index = 0; index < noise_levels; index++){

    double noise = 5;//0.0 + 1 * index;
    Container aux_gp3p(noise, "gp3p");
    Container aux_gpnp(noise, "gpnp");
    Container aux_upnp(noise, "upnp");
    Container aux_nonlinear(noise, "nonlinear");
    Container aux_amm(noise, "amm");
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
				      position, rotation, camOffsets, camRotations                                      , numberPoints, noise, outlierFraction,
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
      
      //add a small perturbation to the rotation
      getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1 );
      adapter.sett(t_perturbed);
      adapter.setR(R_perturbed);
      ObjectiveFunctionInfo * info_container = new GlobalPnPFunctionInfo(adapter);
      double angle =  M_PI / 180;
      double wx = 0; double wy = 0; double wz = 1;
      Eigen::Matrix3d skew_matrix = Eigen::Matrix3d::Zero(3,3);
      skew_matrix(0,1) = -wz; skew_matrix(0,2) = wy;
      skew_matrix(1,0) =  wz; skew_matrix(1,2) = -wx;
      skew_matrix(2,0) = -wy; skew_matrix(2,1) =  wx;
      rotation_t error_rotation = Eigen::Matrix3d::Identity(3,3) +
                              (std::sin(angle) / angle)*skew_matrix +
                              ( ( 1 - std::cos(angle) * std::cos(angle) ) / (angle * angle) ) * skew_matrix * skew_matrix;

      translation_t error_translation = Eigen::Vector3d::Random(3,1);
      error_translation = 0.2 * error_translation / error_translation.norm();
      rotation_t rot = R_perturbed.inverse(); //rotation.inverse();
      translation_t trans = -rot.inverse() * t_perturbed;//position + error_translation;
      //Create solver pointer
      SolverTools * solver_container = NULL;
      solver_container = new SolverToolsNoncentralRelativePose();
      amm solver_object;
      gettimeofday(&tic,0);
      double tol = 1e-6;
      transformation_t amm_solution = solver_object.amm_solver( tol, rot, trans, info_container, solver_container);
      gettimeofday(&toc, 0);
      delete info_container;
      delete solver_container;
     
      double time_amm_solution = TIMETODOUBLE(timeval_minus(toc,tic));
      std::cout << "amm solution: " << std::endl;
      std::cout << "perturbed rotation:      " << std::endl << rot << std::endl;
      std::cout << "perturbed translation:   " << std::endl << trans << std::endl;
      std::cout << "given solution:          " << std::endl << amm_solution << std::endl;
      std::cout << "rotation distance:       " << (rot  - rotation.inverse()).norm() << std::endl;
      std::cout << "translation distance:    " << (trans + rotation.inverse()* position).norm() << std::endl;
      std::cout << "Error rotation (amm):    " << (amm_solution.block<3,3>(0,0) - rotation.inverse()).norm() << std::endl;
      std::cout << "Error translation (amm): " << (amm_solution.block<3,1>(0,3) + rotation.inverse() * position).norm() << std::endl;

      std::cout << "Error rotation (non lin):    " << (nonlinear_transformation.block<3,3>(0,0) - rotation).norm() << std::endl;
      std::cout << "Error translation (non lin): " << (nonlinear_transformation.block<3,1>(0,3) - position).norm() << std::endl;

      std::cout << "Error rotation (gpnp):    " << (gpnp_transformation.block<3,3>(0,0) - rotation).norm() << std::endl;
      std::cout << "Error translation (gpnp): " << (gpnp_transformation.block<3,1>(0,3) - position).norm() << std::endl;
      std::cout << "timings from gp3p algorithm: ";
      std::cout << gp3p_time << std::endl;
      std::cout << "timings from gpnp algorithm: ";
      std::cout << gpnp_time << std::endl;
      std::cout << "timings from upnp algorithm: ";
      std::cout << upnp_time << std::endl;
      std::cout << "timings from nonlinear algorithm: ";
      std::cout << nonlinear_time << std::endl;
      std::cout << "time amm: "                    << time_amm_solution << std::endl;
      int size_upnp = upnp_transformations.size();
      bool flag_non_lin  =  aux_nonlinear.validate_transformation( nonlinear_transformation.block<3,3>(0,0), nonlinear_transformation.block<3,1>(0,3));
      bool flag_gpnp       = aux_gpnp.validate_transformation(gpnp_transformation.block<3,3>(0,0), gpnp_transformation.block<3,1>(0,3));
      bool flag_upnp = aux_upnp.validate_transformation(upnp_transformations[size_upnp - 1].block<3,3>(0,0), upnp_transformations[size_upnp - 1].block<3,1>(0,3));
      bool flag_amm = aux_amm.validate_transformation(amm_solution.block<3,3>(0,0), amm_solution.block<3,1>(0,3));
      if(flag_non_lin == true && flag_gpnp == true && flag_upnp == true &&
	 flag_amm == true){
	aux_nonlinear.add_error_information(rotation, nonlinear_transformation.block<3,3>(0,0), position, nonlinear_transformation.block<3,1>(0,3), nonlinear_time);
	aux_amm.add_error_information(rotation, amm_solution.block<3,3>(0,0),
				      position, amm_solution.block<3,1>(0,3), time_amm_solution);
          aux_gpnp.add_error_information(rotation, gpnp_transformation.block<3,3>(0,0), position, gpnp_transformation.block<3,1>(0,3), gpnp_time);
          aux_upnp.add_error_information(rotation, upnp_transformations[size_upnp - 1].block<3,3>(0,0), position, upnp_transformations[size_upnp - 1].block<3,1>(0,3), upnp_time);
          index_stat++;
      }
      total_realizations++;
    }
    Statistic_info aux(noise, index_stat, total_realizations);
    information_statistics.push_back(aux);
   
    information_gpnp.push_back(aux_gpnp);
    information_upnp.push_back(aux_upnp);
    information_nonlinear.push_back(aux_nonlinear);
    information_amm.push_back(aux_amm);
  }
  std::ofstream datafile;
  datafile.open("data.csv");
  for(unsigned int i = 0; i < information_gpnp.size(); ++i){
    information_gpnp[i].printInfo(datafile);
  }
  for(unsigned int i = 0; i < information_nonlinear.size(); ++i){
    information_nonlinear[i].printInfo(datafile);
  }
  for(unsigned int i = 0; i < information_amm.size(); ++i){
    information_amm[i].printInfo(datafile);
  }
  for(unsigned int i = 0; i < information_upnp.size(); ++i){
    information_upnp[i].printInfo(datafile);
  }
  datafile.close();

  std::ofstream statistical_data;
  statistical_data.open("stat_data.csv");
  for(unsigned int i = 0; i < information_statistics.size(); ++i){
    information_statistics[i].printInfo(statistical_data);
  }
  statistical_data.close();
}
