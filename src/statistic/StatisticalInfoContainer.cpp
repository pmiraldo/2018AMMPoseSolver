#include <iostream>
#include <math.h>
#include <string>
#include <opengv/types.hpp>
#include <opengv/statistic/StatisticalInfoContainer.hpp>
using namespace opengv;

StatisticalInfoContainer::StatisticalInfoContainer(double noise_level_, std::string method_name_, double rotation_error_, double translation_error_, double time_to_run_, std::vector<iterations_info> & information_amm_iterations_) : noise_level(noise_level_), method_name(method_name_), rotation_error(rotation_error_), translation_error(translation_error_), time_to_run(time_to_run_), information_amm_iterations(information_amm_iterations_){
  
}

StatisticalInfoContainer::~StatisticalInfoContainer(){}

void StatisticalInfoContainer::printInfo(std::ostream & stream){

  if (label_error_statistical_info_written == false){
    stream << label_error_statistical_info << std::endl;
    label_error_statistical_info_written = true;
  }
 
    stream << method_name       << ",";
    stream << noise_level       << ",";
    stream << rotation_error    << ",";
    stream << translation_error << ",";
    stream << time_to_run       << std::endl;
 
}

std::string StatisticalInfoContainer::label_error_statistical_info = "method name, noise level, rotation error, translation error, execution time";
bool StatisticalInfoContainer::label_error_statistical_info_written = false;
std::string StatisticalInfoContainer::label_iterations_statistical_info = "method name, noise level, rotation error, translation error, execution time";
bool StatisticalInfoContainer::label_iterations_statistical_info_written = false;
