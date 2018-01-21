#include <opengv/optimization_tools/objective_function_tools/GlobalPnPInfiniteNormFunctionInfo.hpp>
#include <Eigen/Dense>
#include <opengv/Indices.hpp>
#include <iostream>

GlobalPnPInfiniteNormFunctionInfo::GlobalPnPInfiniteNormFunctionInfo(const opengv::absolute_pose::AbsoluteAdapterBase & adapter, const opengv::rotation_t & rotation, const opengv::translation_t &translation){
  Eigen::MatrixXd M   = Eigen::MatrixXd::Zero(3,5);
  Eigen::Matrix3d Vi  = Eigen::MatrixXd::Zero(3,3);
  Eigen::Matrix3d id  = Eigen::MatrixXd::Identity(3,3);
  opengv::Indices indices(adapter.getNumberCorrespondences());
  int total_points = (int) indices.size();
  //std::cout << "Total points: " << total_points << std::endl;
  Eigen::MatrixXd vi     = Eigen::MatrixXd::Zero(3,1);
  Eigen::MatrixXd xi     = Eigen::MatrixXd::Zero(3,1);
  Eigen::MatrixXd ci     = Eigen::MatrixXd::Zero(3,1);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(1,1);

  Eigen::MatrixXd C_all = Eigen::MatrixXd::Zero(3, total_points);
  Eigen::MatrixXd V_all = Eigen::MatrixXd::Zero(3, total_points);
  Eigen::MatrixXd X_all = Eigen::MatrixXd::Zero(3, total_points);
  Eigen::MatrixXd wi    = Eigen::MatrixXd::Zero(3,1);

  Eigen::MatrixXd coefs_gradients = Eigen::MatrixXd::Zero(9,13);
  Eigen::MatrixXd Qi              = Eigen::MatrixXd::Zero(3,3);
  Eigen::MatrixXd Mr              = Eigen::MatrixXd::Zero(3,9);
  for( int i = 0; i < total_points; i++ )
  {
    vi = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
    xi = adapter.getPoint(indices[i]);
    ci = adapter.getCamOffset(indices[i]);
    C_all.block<3,1>(0,i)  = adapter.getCamOffset(indices[i]);
    X_all.block<3,1>(0,i)  = adapter.getPoint(indices[i]);
    V_all.block<3,1>(0,i)  = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
    
    Vi =  vi * vi.transpose()  / ((vi.transpose() * vi)(0,0));
    Qi = (id - Vi).transpose() * (id - Vi);
    M.block<3,3>(0,0) = (id - Vi).transpose() * (id - Vi);
    M.block<3,1>(0,3) = xi;
    M.block<3,1>(0,4) = ci;
    //Coefficients for obj function calculation and the translation
    //gradient are obtained and stored
    coefficients.push_back(M);

    //Initialize list of objective functions
    wi = rotation * xi + translation - ci;
    //std::cout << "wi: " << std::endl << wi << std::endl;
    result = wi.transpose() *  Qi * wi;
    function[i] = result(0,0);

    //Create the matrix used to compute the
    Mr.block<3,3>(0,0) = xi(0,0) * Eigen::MatrixXd::Identity(3,3);
    Mr.block<3,3>(0,3) = xi(1,0) * Eigen::MatrixXd::Identity(3,3);
    Mr.block<3,3>(0,6) = xi(2,0) * Eigen::MatrixXd::Identity(3,3);
    coefs_gradients.block<9,9>(0,0) = Mr.transpose() * Qi * Mr;
    coefs_gradients.block<9,3>(0,9) = 2 * Mr.transpose() * Qi;
    coefs_gradients.block<9,1>(0,12) = -2 * Mr.transpose() * Qi * ci;
    matrix_gradients.push_back(coefs_gradients);
  }

  /*std::cout << "ci"   << std::endl << C_all       << std::endl;
  std::cout << "xi"   << std::endl << X_all       << std::endl;
  std::cout << "vi"   << std::endl << V_all       << std::endl;
  std::cout << "R_: " << std::endl << rotation    << std::endl;
  std::cout << "t_: " << std::endl << translation << std::endl;
  for(int i = 0; i < coefficients.size(); ++i){
    std::cout << "Matrix n: " << i << std::endl << coefficients[i] << std::endl;
  }
  for(int i = 0; i < coefficients.size(); ++i){
    std::cout << "Matrix grad n: " << i << std::endl << matrix_gradients[i] << std::endl;
  }
  for(auto it = function.cbegin(); it != function.cend(); ++it){
    std::cout << it->first << " " << it->second << "\n";
    }*/
 
}

GlobalPnPInfiniteNormFunctionInfo::~GlobalPnPInfiniteNormFunctionInfo(){};

double GlobalPnPInfiniteNormFunctionInfo::objective_function_value(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  Eigen::MatrixXd wi = Eigen::MatrixXd::Zero(3,1);
  Eigen::MatrixXd Qi = Eigen::MatrixXd::Zero(3,3);
  for(std::map<int,double>::iterator it=function.begin(); it != function.end(); ++it){
    wi = rotation * coefficients[it->first].block<3,1>(0,3) + translation - coefficients[it->first].block<3,1>(0,4);
    Qi = coefficients[it->first].block<3,3>(0,0);
    it->second = (wi.transpose() * Qi * wi)(0,0);
  }
   
  using pair_type = decltype(function)::value_type;
  auto pr = std::max_element(
    std::begin(function), std::end(function),
    [] (const pair_type & p1, const pair_type & p2) {
        return p1.second < p2.second;
    }
			     );
  //for(auto it = function.cbegin(); it != function.cend(); ++it){
  //  std::cout << it->first << " " << it->second << "\n";
  //}
  
  min_key = pr->first;
  //std::cout << "MIN KEY: " << min_key << std::endl;
  return (pr->second);
}

opengv::rotation_t GlobalPnPInfiniteNormFunctionInfo::rotation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  objective_function_value(rotation, translation);
  Eigen::MatrixXd Mr  = matrix_gradients[min_key].block<9,9>(0,0);
  Eigen::MatrixXd Mrt = matrix_gradients[min_key].block<9,3>(0,9);
  Eigen::MatrixXd vr  = matrix_gradients[min_key].block<9,1>(0,12);

  const double * p = &rotation(0);
  Map<const Matrix<double,1,9> > r(p, 1, 9);
  Eigen::MatrixXd grad_r = (2 * Mr * r.transpose()) + Mrt * translation + vr;
  double * ptr = &grad_r(0);
  Map<Matrix<double, 3,3> > m(ptr, 3, 3);
  return m;
}

opengv::translation_t GlobalPnPInfiniteNormFunctionInfo::translation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  objective_function_value(rotation, translation);
  Eigen::MatrixXd Qi = coefficients[min_key].block<3,3>(0,0);
  Eigen::MatrixXd xi = coefficients[min_key].block<3,1>(0,3);
  Eigen::MatrixXd ci = coefficients[min_key].block<3,1>(0,4);
  Eigen::MatrixXd grad_t = 2 * Qi * (rotation * xi + translation - ci);
  return grad_t;
}
