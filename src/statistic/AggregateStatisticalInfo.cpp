#include <opengv/statistic/AggregateStatisticalInfo.hpp>
#include <iostream>
#include <math.h>
#include <string>

/*This object concerns the results for a set of experiences for a certain noise level. It is intrinsecally different from the other file which contains information on a single experience. This object merely determines how many trials were needed before reaching the number of trials required*/

AggregateStatisticalInfo::AggregateStatisticalInfo(std::string method_name_, double noise_level_, int successful_realizations_, int total_realizations_) : method_name(method_name_), noise_level(noise_level_), successful_realizations(successful_realizations_), total_realizations(total_realizations_){}

AggregateStatisticalInfo::~AggregateStatisticalInfo(){}

void AggregateStatisticalInfo::printInfo(std::ostream & stream){

  if (label_written == false){
    stream << label << std::endl;
    label_written = true;
  }
  stream << noise_level << ",";
  stream << successful_realizations << ",";
  stream << total_realizations << ",";
  stream << std::endl;
}

std::string AggregateStatisticalInfo::label = "method name, noise level, sucessful, total,";
bool AggregateStatisticalInfo::label_written = false;
