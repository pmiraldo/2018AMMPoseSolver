#ifndef AMM_H
#define AMM_H

#include <opengv/optimization_tools/solver_tools/SolverTools.hpp>
#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>
#include <opengv/statistic/iterations_info.hpp>
#include <opengv/types.hpp>
#include <vector>



class amm{
public:
  amm();
  ~amm();
  opengv::transformation_t amm_solver(double & tol,
				      const opengv::rotation_t & initial_state,
				      const opengv::translation_t & initial_translation,
				      ObjectiveFunctionInfo * objective_function_container,
				      SolverTools * solver_container, double & step, std::vector<iterations_info> & iter_stat_list);
  std::vector<iterations_info> get_iterations_result(){return list;}
private:
  std::vector<iterations_info> list;
}; 
#endif
