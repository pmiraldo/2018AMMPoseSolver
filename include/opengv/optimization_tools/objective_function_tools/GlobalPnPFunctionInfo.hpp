#ifndef GLOBALPNPFUNCTIONINFO_H
#define GLOBALPNPFUNCTIONINFO_H

#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>
#include <opengv/absolute_pose/AbsoluteAdapterBase.hpp>
#include <opengv/types.hpp>

class GlobalPNPFunctionInfo : public ObjectiveFunctionInfo {
public:
  GlobalPNPFunctionInfo(const opengv::absolute_pose::AbsoluteAdapterBase & adapter );
  ~GlobalPNPFunctionInfo();

  /*  double GlobalPNPFunctionInfo::objective_function_value(const opengv::rotation_t & rotation,
                                          const opengv::translation_t & translation);
  opengv::rotation_t GlobalPNPFunctionInfo::rotation_gradient(const opengv::rotation_t & rotation,
							      const opengv::translation_t & translation);
  opengv::translation_t GlobalPNPFunctionInfo::translation_gradient(const opengv::rotation_t & rotation,
  const opengv::translation_t & translation);*/

private:
  Eigen::Matrix3d Mt;
  Eigen::MatrixXd Mrt;
  Eigen::MatrixXd Mr; 

  Eigen::Vector3d vt; 
  Eigen::VectorXd vr;

};

#endif
