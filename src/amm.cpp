#include <opengv/amm.hpp>
#include <opengv/types.hpp>
#include <opengv/optimization_tools/solver_tools/SolverTools.hpp>
#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>
#include <opengv/statistic/iterations_info.hpp>
#include <iostream>
#include <vector>


amm::amm(){}

amm::~amm(){}

opengv::transformation_t amm::amm_solver(double & tol,
                                      const opengv::rotation_t & initial_state,
                                      const opengv::translation_t & initial_translation,
                                      ObjectiveFunctionInfo * objective_function_container,
					 SolverTools * solver_container, double & step, std::vector<iterations_info> & iter_stat_list){
  double error = 1;
  int max = 100;
  int iteration = 0;
  opengv::rotation_t state = initial_state;
  opengv::translation_t translation = initial_translation;
  //opengv::rotation_t previous_state = initial_state;
  double tol_solvers = 1e-6;
  //std::cout << "Beginning of the amm" << std::endl;
  int translation_iterations = 0;
  int rotation_iterations = 0;
  iter_stat_list.clear();
  iterations_info aux;
  while (error > tol && iteration < max){
    //previous_state = state;
    //minimize rotation state
    //std::cout << "Iteration: " << iteration << std::endl;
    double f_previous_value = objective_function_container->objective_function_value(state, translation);
    translation = solver_container->translation_solver(state, translation, tol_solvers, objective_function_container, step, translation_iterations);
    /*std::cout << "New rotation: "              << std::endl << state       << std::endl;
      std::cout << "The translation: "           << std::endl << translation << std::endl;*/
    state = solver_container->rotation_solver(state, translation, tol_solvers, objective_function_container, rotation_iterations);
      /*std::cout << "New translation: "           << std::endl << translation << std::endl;
	std::cout << "End of iteration: " << std::endl;*/
   double f_current_value = objective_function_container->objective_function_value(state, translation);
   error = std::abs(f_previous_value - f_current_value);
   iteration++;
  
   aux.iterations_rotation = rotation_iterations;
   aux.iterations_translation = translation_iterations;
   iter_stat_list.push_back(aux);
  }
 
  opengv::transformation_t solution;
  solution.block<3,3>(0,0) = state;
  solution.block<3,1>(0,3) = translation;
  return solution;
}


