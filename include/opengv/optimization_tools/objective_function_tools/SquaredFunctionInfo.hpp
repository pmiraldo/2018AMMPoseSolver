#ifndef SQUAREDFUNCTIONINFO_H
#define SQUAREDFUNCTIONINFO_H

#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>
#include <Eigen/Dense>
#include <opengv/types.hpp>
#include <opengv/relative_pose/RelativeAdapterBase.hpp>

class SquaredFunctionInfo : public ObjectiveFunctionInfo {
  Eigen::MatrixXd M;

 public:
  SquaredFunctionInfo(const opengv::relative_pose::RelativeAdapterBase & adapter );
  ~SquaredFunctionInfo();

  Eigen::MatrixXd get_M();
  double objective_function_value(const opengv::rotation_t & rotation, const opengv::translation_t & translation);
  opengv::rotation_t rotation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation);
  opengv::translation_t translation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation);
};
		    
#endif
