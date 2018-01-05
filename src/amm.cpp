#include <opengv/amm.hpp>
#include <opengv/types.hpp>
#include <opengv/optimization_tools/solver_tools/SolverTools.hpp>
#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>

amm::amm(){}

amm::~amm(){}

opengv::transformation_t amm::amm_solver(double & tol,
                                      const opengv::rotation_t & initial_state,
                                      const opengv::translation_t & initial_translation,
                                      ObjectiveFunctionInfo * objective_function_container,
				      SolverTools * solver_container){
  double error = 1;
  int max = 100;
  int iteration = 0;
  opengv::rotation_t state = initial_state;
  opengv::translation_t translation = initial_translation;
  opengv::rotation_t previous_state = initial_state;
  double tol_solvers = 1e-6;
  while (error > tol && iteration < max){
    previous_state = state;
    //minimize rotation state
    /*std::cout << "**************************************************"      << std::endl;
    std::cout << "Iteration: "                              << iteration   << std::endl;
    std::cout << "Rotation at the beginning: " << std::endl << state       << std::endl;*/
    state = solver_container->rotation_solver(state, translation, tol_solvers, objective_function_container);

    /*std::cout << "New rotation: "              << std::endl << state       << std::endl;
    std::cout << "The translation: "           << std::endl << translation << std::endl;*/
    translation = solver_container->translation_solver(state, translation, tol_solvers, objective_function_container);
    /*std::cout << "New translation: "           << std::endl << translation << std::endl;
      std::cout << std::endl << std::endl;*/
    error = (previous_state - state).norm();
    iteration++;
  }
 
  opengv::transformation_t solution;
  solution.block<3,3>(0,0) = state;
  solution.block<3,1>(0,3) = translation;
  return solution;
}


