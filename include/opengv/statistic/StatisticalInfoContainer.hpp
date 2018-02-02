#ifndef STATISTICALINFOCONTAINER_HPP
#define STATISTICALINFOCONTAINER_HPP

#include <vector>
#include <opengv/types.hpp>
#include <opengv/statistic/iterations_info.hpp>
#include <string>
#include <iostream>
#include <fstream>

class StatisticalInfoContainer{
 public:
  StatisticalInfoContainer(double noise_level_, std::string method_name_, double rotation_error_, double translation_error_, double time_to_run_, std::vector<iterations_info> & information_amm_iterations_);
  ~StatisticalInfoContainer();
  
  void printInfo(std::ostream & stream);

 private:
  std::string method_name;
  double noise_level;
  double rotation_error;
  double translation_error;
  double time_to_run;
  std::vector<iterations_info> information_amm_iterations;

  //static members to control writing the label
  static std::string label_error_statistical_info;
  static bool label_error_statistical_info_written;
  
  static std::string label_iterations_statistical_info;
  static bool label_iterations_statistical_info_written;


  double calc_rotation_error(const opengv::rotation_t & original_rotation,
			     const opengv::rotation_t & obtained_rotation);

  double calc_translation_error(const opengv::translation_t & original_translation,
				const opengv::translation_t & obtained_translation);

};

#endif
