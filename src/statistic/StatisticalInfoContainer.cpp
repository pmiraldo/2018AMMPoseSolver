#include <iostream>
#include <math.h>
#include <string>
#include <opengv/types.hpp>
#include <opengv/statistic/StatisticalInfoContainer.hpp>
using namespace opengv;

StatisticalInfoContainer::StatisticalInfoContainer(double noise_level_, std::string method_name_, double rotation_error_, double translation_error_, double time_to_run_, std::vector<iterations_info> & information_amm_iterations_) : noise_level(noise_level_), method_name(method_name_), rotation_error(rotation_error_), translation_error(translation_error_), time_to_run(time_to_run_), information_amm_iterations(information_amm_iterations_){
  
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
