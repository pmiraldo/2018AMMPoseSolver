#include "container.h"
#include <iostream>
#include <math.h>
#include <string>

using namespace opengv;

Container::Container(double noise_level_, std::string method_name_){
  noise_level = noise_level_;
  method_name = method_name_;
}

Container::~Container(){}

bool Container::validate_transformation(const opengv::rotation_t    & obtained_rotation,
			     const opengv::translation_t & obtained_translation){
  
  bool flag00 = std::isnan(obtained_rotation(0,0));
  bool flag01 = std::isnan(obtained_rotation(0,1));
  bool flag02 = std::isnan(obtained_rotation(0,2));
  bool flag10 = std::isnan(obtained_rotation(1,0));
  bool flag11 = std::isnan(obtained_rotation(1,1));
  bool flag12 = std::isnan(obtained_rotation(1,2));
  bool flag20 = std::isnan(obtained_rotation(2,0));
  bool flag21 = std::isnan(obtained_rotation(2,1));
  bool flag22 = std::isnan(obtained_rotation(2,2));

  bool flagt01 = std::isnan(obtained_translation(0,0));
  bool flagt02 = std::isnan(obtained_translation(0,1));
  bool flagt03 = std::isnan(obtained_translation(0,2));
  if (flag00  == true || flag01  == true || flag02  == true ||
      flag10  == true || flag11  == true || flag12  == true ||
      flag20  == true || flag21  == true || flag22  == true ||
      flagt01 == true || flagt02 == true || flagt03 == true){
    return false;
  }
  return true;
}
void Container::add_error_information(const rotation_t    & original_rotation   , const rotation_t    & obtained_rotation,
				      const translation_t & original_translation, const translation_t & obtained_translation,
				      const double & time){
  //check if the obtained rotation is a valid one:
  bool flag = validate_transformation(obtained_rotation, obtained_translation);
  if (flag == false){
    return;
  }
  /* If it passed this stage is because it is a valid matrix. It is assumed that after this, every information is valid	\
     .*/
  info_node aux;
  aux.error_rotation    = calc_rotation_error(original_rotation, obtained_rotation);
  aux.error_translation = calc_translation_error(original_translation, obtained_translation);
  aux.execution_time    = time;
  statistic_info.push_back(aux);
}

void Container::printInfo(std::ostream & stream){

  if (label_written == false){
    stream << label << std::endl;
    label_written = true;
  }
  for(unsigned int i = 0; i < statistic_info.size(); ++i){
    stream << method_name << ",";
    stream << noise_level << ",";
    stream << statistic_info[i].error_rotation    << ",";
    stream << statistic_info[i].error_translation << ",";
    stream << statistic_info[i].execution_time    << std::endl;
  }
}

double Container::calc_rotation_error(const rotation_t & original_rotation, const rotation_t & obtained_rotation){
  return( (original_rotation - obtained_rotation).norm() );
}

double Container::calc_translation_error(const translation_t & original_translation, const translation_t & obtained_translation){
  return( (original_translation - obtained_translation).norm() );
}

std::string Container::label = "method name, noise level, rotation error, translation error, execution time";
bool Container::label_written = false;
