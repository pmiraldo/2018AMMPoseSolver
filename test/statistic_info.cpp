#include "statistic_info.h"
#include <iostream>
#include <math.h>
#include <string>

using namespace opengv;

Statistic_info::Statistic_info(double noise_level_, int success_, int total_){
  noise_level = noise_level_;
  successful_realizations = success_;
  total_realizations = total_;
}

Statistic_info::~Statistic_info(){}


void Statistic_info::printInfo(std::ostream & stream){

  if (label_written == false){
    stream << label << std::endl;
    label_written = true;
  }
  stream << noise_level << ",";
  stream << successful_realizations << ",";
  stream << total_realizations << ",";
  stream << std::endl;
}

std::string Statistic_info::label = "noise level, sucessful, total,";
bool Statistic_info::label_written = false;
