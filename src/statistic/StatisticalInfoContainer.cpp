#include <iostream>
#include <math.h>
#include <string>
#include <opengv/types.hpp>
#include <opengv/statistic/StatisticalInfoContainer.hpp>
using namespace opengv;

StatisticalInfoContainer::StatisticalInfoContainer(double noise_level_, std::string method_name_, const opengv::transformation_t & transformation_method, const opengv::transformation_t & transformation_reference, double time_to_run_, std::vector<iterations_info> & information_amm_iterations_) : method_name(method_name_), noise_level(noise_level_), time_to_run(time_to_run_), information_amm_iterations(information_amm_iterations_)
{
  rotation_error = (transformation_method.block<3,3>(0,0) - transformation_reference.block<3,3>(0,0)).norm();
  translation_error = (transformation_method.block<3,1>(0,3) - transformation_reference.block<3,1>(0,3)).norm();
}

StatisticalInfoContainer::~StatisticalInfoContainer(){}

void StatisticalInfoContainer::printInfo(std::ostream & error_information, std::ostream & iteration_information, const bool & is_amm){

  if (label_error_statistical_info_written == false){
    error_information << label_error_statistical_info << std::endl;
    label_error_statistical_info_written = true;
  }
 
  error_information << method_name       << ",";
  error_information << noise_level       << ",";
  error_information << rotation_error    << ",";
  error_information << translation_error << ",";
  error_information << time_to_run       << std::endl;
  if (is_amm == true){
    if (label_iterations_statistical_info_written == false){
      iteration_information << label_iterations_statistical_info << std::endl;
      label_iterations_statistical_info_written = true;
    }
    iteration_information << method_name       << ",";
    iteration_information << noise_level       << ",";
    iteration_information << "rotation: ";
    for(int i = 0; i < information_amm_iterations.size(); ++i){
      
      iteration_information << information_amm_iterations[i].iterations_rotation;
      if( i != information_amm_iterations.size() - 1){
	iteration_information <<  " ,";
      }
    }
    iteration_information << std::endl;
    iteration_information << method_name       << ",";
    iteration_information << noise_level       << ",";
    iteration_information << "translation: ";
    for(int i = 0; i < information_amm_iterations.size(); ++i){
      iteration_information << information_amm_iterations[i].iterations_translation;
      if( i != information_amm_iterations.size() - 1){
	iteration_information <<  " ,";
      }
    }
    iteration_information << std::endl;
  }
}

std::string StatisticalInfoContainer::label_error_statistical_info = "method name, noise level, rotation error, translation error, execution time";
bool StatisticalInfoContainer::label_error_statistical_info_written = false;
std::string StatisticalInfoContainer::label_iterations_statistical_info = "method name, noise level, iteration kind";
bool StatisticalInfoContainer::label_iterations_statistical_info_written = false;
