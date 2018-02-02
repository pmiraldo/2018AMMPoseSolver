#ifndef AGGREGATESTATISTICALINFO_HPP
#define AGGREGATESTATISTICALINFO_HPP

#include <vector>
#include <opengv/types.hpp>
#include <string>
#include <iostream>
#include <fstream>

class AggregateStatisticalInfo{
 public:
  AggregateStatisticalInfo(std::string method_name_, double noise_level_, int successful_realizations_, int total_realizations_);
  ~AggregateStatisticalInfo();
  
  void printInfo(std::ostream & stream);

 private:
  double noise_level;
  int successful_realizations;
  int total_realizations;
  std::string method_name;
  //static members to control writting the label
  static std::string label;
  static bool label_written;
};
#endif
