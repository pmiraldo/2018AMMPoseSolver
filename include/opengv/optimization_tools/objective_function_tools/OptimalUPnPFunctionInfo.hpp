#ifndef OPTIMALUPNPFUNCTIONINFO_H
#define OPTIMALUPNPFUNCTIONINFO_H

#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>
#include <opengv/absolute_pose/AbsoluteAdapterBase.hpp>
#include <opengv/types.hpp>

class OptimalUPnPFunctionInfo : public ObjectiveFunctionInfo {
public:
  OptimalUPnPFunctionInfo(const opengv::absolute_pose::AbsoluteAdapterBase & adapter, opengv::rotation_t rot, opengv::translation_t trans);
  ~OptimalUPnPFunctionInfo();

  double objective_function_value(const opengv::rotation_t & rotation,
                                          const opengv::translation_t & translation);
  opengv::rotation_t rotation_gradient(const opengv::rotation_t & rotation,
							      const opengv::translation_t & translation);
  opengv::translation_t translation_gradient(const opengv::rotation_t & rotation,
  const opengv::translation_t & translation);

private:
  Eigen::MatrixXd Mr;
  Eigen::MatrixXd vr;
  Eigen::MatrixXd Mrt;
  Eigen::MatrixXd vt;
  int n;
};

#endif
