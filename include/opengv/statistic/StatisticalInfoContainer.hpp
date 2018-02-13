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
  StatisticalInfoContainer(int n_trials_, int noise_levels_, int n_iterations_, std::vector<std::string> methods_list_);

  ~StatisticalInfoContainer();

  void append_info(int method_index, int index_n, double noise, const opengv::transformation_t & method_result, const opengv::transformation_t & gt_transformation, const double & execution_time, std::vector<iterations_info> number_of_interations);
  
  void print_info(std::string filename);

 private:
  std::vector<std::string> methods_list;
  Eigen::MatrixXd trial_info_table;
  int block_size;
  int n_iterations;
};

#endif
