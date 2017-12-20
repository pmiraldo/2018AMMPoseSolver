#ifndef CONTAINER_H
#define CONTAINER_H

#include <vector>
#include <opengv/types.hpp>
#include <string>
#include <iostream>
#include <fstream>

struct info_node{
  double error_rotation;
  double error_translation;
  double execution_time;
};

class Container{
 public:
  Container(double noise_level_, std::string method_name_);
  ~Container();
  bool validate_transformation(const opengv::rotation_t    & obtained_rotation,
			       const opengv::translation_t & obtained_translation);
  
  void add_error_information(const opengv::rotation_t    & original_rotation,
			     const opengv::rotation_t    & obtained_rotation,
			     const opengv::translation_t & original_translation,
			     const opengv::translation_t & obtained_translation,
			     const double & time);
  double get_noise_level(){return noise_level;}

  void printInfo(std::ostream & stream);

 private:
  std::string method_name;
  double noise_level;
  std::vector<info_node> statistic_info;
  //static members to control writing the label
  static std::string label;
  static bool label_written;

  double calc_rotation_error(const opengv::rotation_t & original_rotation,
			     const opengv::rotation_t & obtained_rotation);

  double calc_translation_error(const opengv::translation_t & original_translation,
				const opengv::translation_t & obtained_translation);

};

#endif
