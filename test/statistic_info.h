#ifndef STATISTIC_INFO_H
#define STATISTIC_INFO_H

#include <vector>
#include <opengv/types.hpp>
#include <string>
#include <iostream>
#include <fstream>

class Statistic_info{
 public:
  Statistic_info(double noise_level_, int success_, int total_);
  ~Statistic_info();
  
  void printInfo(std::ostream & stream);

 private:
  double noise_level;
  int successful_realizations;
  int total_realizations;
  //static members to control writting the label
  static std::string label;
  static bool label_written;
};
#endif
