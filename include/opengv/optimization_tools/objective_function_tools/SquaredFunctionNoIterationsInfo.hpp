#ifndef SQUAREDFUNCTIONNOITERATIONSINFO_H
#define SQUAREDFUNCTIONNOITERATIONSINFO_H

#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>
#include <Eigen/Dense>
#include <opengv/types.hpp>
#include <opengv/relative_pose/RelativeAdapterBase.hpp>

class SquaredFunctionNoIterationsInfo : public ObjectiveFunctionInfo {
  Eigen::MatrixXd M;
  Eigen::MatrixXd v;
 public:
  SquaredFunctionNoIterationsInfo(const opengv::relative_pose::RelativeAdapterBase & adapter );
  ~SquaredFunctionNoIterationsInfo();

  Eigen::MatrixXd get_M();
  double objective_function_value(const opengv::rotation_t & rotation, const opengv::translation_t & translation);
  opengv::rotation_t rotation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation);
  opengv::translation_t translation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation);
};
		    
#endif
