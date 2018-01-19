#ifndef GLOBALPNPINFINITENORMFUNCTIONINFO_H
#define GLOBALPNPINFINITENORMFUNCTIONINFO_H

#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>
#include <opengv/absolute_pose/AbsoluteAdapterBase.hpp>
#include <opengv/types.hpp>
#include <vector>
#include <map>

class GlobalPnPInfiniteNormFunctionInfo : public ObjectiveFunctionInfo {
public:
  GlobalPnPInfiniteNormFunctionInfo(const opengv::absolute_pose::AbsoluteAdapterBase & adapter, const opengv::rotation_t & rotation, const opengv::translation_t & translation);
  ~GlobalPnPInfiniteNormFunctionInfo();

  double objective_function_value(const opengv::rotation_t & rotation,
				  const opengv::translation_t & translation);
  opengv::rotation_t rotation_gradient(const opengv::rotation_t & rotation,
				       const opengv::translation_t & translation);
  opengv::translation_t translation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation);

private:
  std::vector<Eigen::MatrixXd> coefficients;
  std::map<int, double> function;
  int min_key = 0;
  std::vector<Eigen::MatrixXd> matrix_gradients;
};

#endif

