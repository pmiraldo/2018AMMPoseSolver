#ifndef AMM_H
#define AMM_H

#include <opengv/optimization_tools/solver_tools/SolverTools.hpp>
#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>
#include <opengv/types.hpp>


class amm{
public:
  amm();
  ~amm();
  opengv::transformation_t amm_solver(double & tol,
				      const opengv::rotation_t & initial_state,
				      const opengv::translation_t & initial_translation,
				      ObjectiveFunctionInfo * objective_function_container,
				      SolverTools * solver_container);
 
}; 
#endif
