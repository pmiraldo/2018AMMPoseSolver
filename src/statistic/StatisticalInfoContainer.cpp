#include <iostream>
#include <math.h>
#include <string>
#include <opengv/types.hpp>
#include <opengv/statistic/StatisticalInfoContainer.hpp>
using namespace opengv;

StatisticalInfoContainer::StatisticalInfoContainer(int n_trials_, int noise_levels_, int n_iterations_, std::vector<std::string> methods_list_) :  n_iterations(n_iterations_),  methods_list(methods_list_) {

  int total_rows = n_trials_ * noise_levels_ * methods_list.size();
  int total_cols = 5 + 2 * n_iterations; //5: method name, noise level, rotatior error, translation error, execution time + 20 rotation solver iterations + 20 translation solver iterations
  block_size = n_trials_ * noise_levels_;
  trial_info_table = Eigen::MatrixXd::Zero(total_rows, total_cols);
  std::cout << "The table was created: " << std::endl;
  std::cout << "cols: " << trial_info_table.cols() << std::endl;
  std::cout << "rows: " << trial_info_table.rows() << std::endl;
}
  

void StatisticalInfoContainer::append_info(int method_index, int index_n, double noise, const opengv::transformation_t & method_result, const opengv::transformation_t & gt_transformation, const double & execution_time, std::vector<iterations_info> number_of_iterations){

  int row_to_update = block_size * method_index + index_n;
  trial_info_table(row_to_update, 0) = method_index;
  trial_info_table(row_to_update, 1) = noise;
  trial_info_table(row_to_update, 2) = (gt_transformation.block<3,3>(0,0) - method_result.block<3,3>(0,0) ).norm();
  trial_info_table(row_to_update, 3) = (gt_transformation.block<3,1>(0,3) - method_result.block<3,1>(0,3) ).norm();
  trial_info_table(row_to_update, 4) = execution_time;
  if(number_of_iterations.size() > 0){
    int max_index = (n_iterations > number_of_iterations.size()) ? (number_of_iterations.size()) : (n_iterations);
    int offset = 5;
    for(int i = 0; i < max_index; ++i){
      trial_info_table(row_to_update, offset + i) = number_of_iterations[i].iterations_rotation;
      trial_info_table(row_to_update, offset + i + n_iterations) = number_of_iterations[i].iterations_translation;
    }
    }
  std::cout << "Row corresponding to method: " << methods_list[method_index] << " " << row_to_update << std::endl;
}
					   
StatisticalInfoContainer::~StatisticalInfoContainer(){
}

void StatisticalInfoContainer::print_info(std::string filename){
  std::ofstream  statistical_information;
  statistical_information.open(filename);

  std::string label = "method_name,noise_levels,rotation_error,translation_error,execution_time";
  for(int i = 0; i < n_iterations; ++i){
    label.append(",exp_rot_" + std::to_string(i));
  }
  for(int i = 0; i < n_iterations; ++i){
    label.append(",exp_trans_" + std::to_string(i));
  }
  
  statistical_information << label << std::endl;
  std::cout << "Label: " << label << std::endl;
  for(int i = 0; i < trial_info_table.rows(); ++i){
    statistical_information << methods_list[(int) trial_info_table(i,0)];
    for(int j = 1; j < trial_info_table.cols(); ++j){
      statistical_information << "," << trial_info_table(i,j);
    }
    statistical_information << std::endl;
  }
  statistical_information.close();
}

