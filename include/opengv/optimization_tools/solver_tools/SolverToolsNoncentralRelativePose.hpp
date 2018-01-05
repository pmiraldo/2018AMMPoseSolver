#include <opengv/optimization_tools/solver_tools/SolverTools.hpp>
#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>
#include <Eigen/Dense>
#include <opengv/types.hpp>

class SolverToolsNoncentralRelativePose : public SolverTools {
public:
  Eigen::Matrix3d exp_R( Eigen::Matrix3d & X );

  opengv::rotation_t rotation_solver(opengv::rotation_t & state_rotation, const opengv::translation_t & translation,
                                     double &tol, ObjectiveFunctionInfo * info_function);

  opengv::translation_t translation_solver(const opengv::rotation_t & rotation,
                                           opengv::translation_t & translation, double &tol, ObjectiveFunctionInfo * info_function);

};
