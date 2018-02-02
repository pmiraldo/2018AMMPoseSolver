#include <opengv/optimization_tools/objective_function_tools/ObjectiveFunctionInfo.hpp>
#include <opengv/optimization_tools/solver_tools/SolverToolsNoncentralRelativePose.hpp>
#include <Eigen/Dense>
#include <opengv/types.hpp>
#include <iostream>


Eigen::Matrix3d SolverToolsNoncentralRelativePose::exp_R( Eigen::Matrix3d & X ){
  double phi = X.norm()/std::sqrt(2);
  Eigen::Matrix3d X1 = X/phi;
  Eigen::Matrix3d I = Eigen::Matrix< double, 3, 3 >::Identity();
  I = I + std::sin(phi)*X1 + ( 1 - std::cos(phi) )*X1*X1;
  return I;
}


opengv::rotation_t SolverToolsNoncentralRelativePose::rotation_solver(opengv::rotation_t & state_rotation, const opengv::translation_t & translation,
								      double &tol, ObjectiveFunctionInfo * info_function, int & k){
  double g = 1.0;
  double erro = 1.0;
  k = 0;
  opengv::rotation_t X = state_rotation;
  opengv::rotation_t previous_X = Eigen::Matrix3d::Zero(3,3);
  Eigen::Matrix3d Z  = Eigen::Matrix3d::Zero(3,3);
  double zz = 0;
  Eigen::Matrix3d DX = Eigen::Matrix3d::Zero(3,3);
  Eigen::Matrix3d Pt = Eigen::Matrix3d::Zero(3,3);
  Eigen::Matrix3d P  = Eigen::Matrix3d::Zero(3,3);
  Eigen::Matrix3d Q  = Eigen::Matrix3d::Zero(3,3);
  Eigen::Matrix3d Qt = Eigen::Matrix3d::Zero(3,3);
  Eigen::Matrix3d reference = Eigen::Matrix3d::Identity(3,3);
  //std::cout << "Inside rotation solver to check the values: " << std::endl;
  //std::cout << "Before process starts the rotation matrix is: " << X << std::endl;
  
  while( erro > tol && k < 1e5 )
    {

      //std::cout << "state at the beginning: " << std::endl << X << std::endl << std::endl;
      DX = info_function->rotation_gradient(X, translation);
      //std::cout << "Calculated gradient: " << std::endl << DX << std::endl << std::endl;
      Z   = DX*X.transpose() - X*DX.transpose();
      //std::cout << "Calculated riemaniann gradient: " << std::endl << Z << std::endl << std::endl;
      zz  = 0.5*( Z*Z.transpose() ).trace();
      //std::cout << "Coefficient zz : " << zz << std::endl << std::endl;
      Pt  = -g*Z;
      //std::cout << "Matrix Pt: " << std::endl << Pt << std::endl << std::endl;
      P   = exp_R( Pt );
      //std::cout << "The rotation matrix P: " << std::endl << P << std::endl << std::endl;
      Q   = P*P; // this seems strange
      //std::cout << "Matrix Q: " << std::endl << Q << std::endl << std::endl;
      Qt  = Q*X;
      //std::cout << "Matrix Qt: " << std::endl << Qt << std::endl << std::endl;
      //std::cout << "************************************************************" << std::endl << std::endl;
      //std::cout << "1st CYCLE" << std::endl << std::endl;
      
      //while( ( objective_function( M, X, translation ) - objective_function(M, Qt, translation ) ) >= g*zz  )
      while( ( info_function->objective_function_value(X, translation ) - info_function->objective_function_value(Qt, translation ) ) >= g*zz  )
        {
          g   = 2*g;
          //In order to prevent NAN's the following restriction is added
          if(g < 64){
            //std::cout << "g: " << g << std::endl;
            P   = Q;
            //std::cout << std::endl << "New P: " << std::endl << P << std::endl;
            Q   = P*P; // this seems strange
            //std::cout << std::endl << "New Q: " << std::endl << Q << std::endl;
            Qt  = Q*X;
            //std::cout << std::endl << "New Qt: " << std::endl << Qt << std::endl;
            //std::cout << "Current value for obj function: " << objective_function(M, X, translation) << std::endl;
            //std::cout << "Current value for new obj function: " << objective_function(M, Qt, translation ) << std::endl;
          }
          else{
            break;
          }
        }
      //std::cout << "End of 1st cycle" << std::endl;
      //std::cout << "******************************************************************" << std::endl;
      //std::cout << "**********************************************************************" << std::endl;
      //std::cout << "Enters 2nd cycle: " << std::endl;
      Qt = P * X;
      //while( ( objective_function( M, X, translation ) - objective_function(M, Qt, translation ) ) < 0.5*g*zz)
      while( ( info_function->objective_function_value( X, translation ) - info_function->objective_function_value( Qt, translation ) ) < 0.5*g*zz)
        {

          //   if ( f_obj( M, N, X, beta) - f_obj( M, N, Qt, beta ) < tol )
          //     break;

          g  = 0.5*g;
          //std::cout << "New g: " << std::endl << g << std::endl << std::endl;
          Pt = -g*Z;
          //std::cout << "New Pt: " << std::endl << Pt << std::endl << std::endl;
          P  = exp_R( Pt );
          //In order to prevent NAN's
          if( (P - reference).norm() < 1e-6){
            break;
          }
          //std::cout << "New P : " << std::endl << P << std::endl << std::endl;
          Qt = P*X;
          //std::cout << "New Qt: " << std::endl << std::endl << Qt << std::endl;
        }
      previous_X = X;
      X    = P*X;
      erro = ( X - previous_X ).norm();
      /*std::cout << "inside rotation solver: " << std::endl;
      std::cout << "\nIteration: " << k << std::endl;
      std::cout << "Rotation: "    << std::endl << X << std::endl;
      std::cout << "Euclidean grad: " << std::endl << DX << std::endl;
      std::cout << "Function value: " << info_function->objective_function_value(X, translation) << std::endl;*/
      k++;
    }
  return X;
}



opengv::translation_t SolverToolsNoncentralRelativePose::translation_solver(const opengv::rotation_t & rotation, opengv::translation_t & translation, double &tol, ObjectiveFunctionInfo * info_function, double & step, int & k){

 
  double error = 1;
  k = 0;
  //std::cout << "Translation gradient: " << std::endl;
  opengv::translation_t state = translation;
  opengv::translation_t grad = info_function->translation_gradient(rotation, state);
  opengv::translation_t new_state = state - step * grad;
  //std::cout << "Beginning of the translation solver: " << std::endl;
  //std::cout << "The rotation used is: " << std::endl << rotation << std::endl;
  while (error > tol ){
    new_state = state - step * grad;
    double f_obj_current = info_function->objective_function_value(rotation, state);
    double f_obj_next = info_function->objective_function_value(rotation, new_state);

    /*std::cout << "current state: " << std::endl     << state      << std::endl;
    std::cout << "new state: "     << std::endl     << new_state  << std::endl;
    std::cout << "f(state): "      << f_obj_current << std::endl;
    std::cout << "f(new_state): "  << f_obj_next    << std::endl;
    std::cout << "gradient:     "  << std::endl     << grad       << std::endl;
    std::cout << std::endl         << std::endl     << std::endl;*/
    if(f_obj_next > f_obj_current){
      break;
    }
    error = std::abs(f_obj_current - f_obj_next);
    grad = info_function->translation_gradient(rotation, state);
    
    state = new_state;
    k++;
  }
  return state;
}

