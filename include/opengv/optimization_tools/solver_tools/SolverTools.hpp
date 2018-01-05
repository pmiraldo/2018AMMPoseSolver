#ifndef SOLVERTOOLS_H
#define SOLVERTOOLS_H

#include <opengv/types.hpp>
#include <Eigen/Dense>
#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>

class ObjectiveFunctionInfo;

class SolverTools{
 public:

  virtual Eigen::Matrix3d exp_R( Eigen::Matrix3d & X ) = 0;
 
  virtual opengv::rotation_t rotation_solver(opengv::rotation_t & state_rotation, const opengv::translation_t & translation,
					     double &tol, ObjectiveFunctionInfo * info_function) = 0;

  virtual opengv::translation_t translation_solver(const opengv::rotation_t & rotation,
						   opengv::translation_t & translation, double &tol, ObjectiveFunctionInfo * info_function) = 0;

};

#endif
