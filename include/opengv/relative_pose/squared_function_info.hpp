#ifndef SQUARED_FUNCTION_INFO_H
#define SQUARED_FUNCTION_INFO_H

#include <opengv/relative_pose/objective_function_info.hpp>
#include <Eigen/Dense>
#include <opengv/relative_pose/RelativeAdapterBase.hpp>


class squared_function_info : public objective_function_info {
  Eigen::MatrixXd M;

 public:
  squared_function_info(const opengv::relative_pose::RelativeAdapterBase & adapter );
  ~squared_function_info();

  Eigen::MatrixXd get_M();
  double objective_function_value(const opengv::rotation_t & rotation, const opengv::translation_t & translation);
  opengv::rotation_t rotation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation);
  opengv::translation_t translation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation);
};
		    
#endif
