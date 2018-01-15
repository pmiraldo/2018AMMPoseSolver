#include <opengv/optimization_tools/objective_function_tools/OptimalUPnPFunctionInfo.hpp>
#include <Eigen/Dense>
#include <opengv/Indices.hpp>
#include <iostream>

OptimalUPnPFunctionInfo::OptimalUPnPFunctionInfo(const opengv::absolute_pose::AbsoluteAdapterBase & adapter, opengv::rotation_t rot, opengv::translation_t trans){
  //Initialize class members
  Mr  = Eigen::MatrixXd::Zero(9,9);
  vr  = Eigen::MatrixXd::Zero(9,1);;
  Mrt = Eigen::MatrixXd::Zero(9,3);
  vt  = Eigen::MatrixXd::Zero(3,1);
  Eigen::MatrixXd constant = Eigen::MatrixXd::Zero(1,1); //Used for debugging
  n   = 0;
  opengv::Indices indices(adapter.getNumberCorrespondences());
  int total_points = (int) indices.size();
  n = total_points;
  
  //Used to store all the information needed
  Eigen::MatrixXd C_all = Eigen::MatrixXd::Zero(3, total_points);
  Eigen::MatrixXd V_all = Eigen::MatrixXd::Zero(3, total_points);
  Eigen::MatrixXd X_all = Eigen::MatrixXd::Zero(3, total_points);
  //**************************
  Eigen::Matrix3d id = Eigen::Matrix3d::Identity(3,3);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3 * n, n + 3);
  // This time all the information is needed beforehand
  for( int i = 0; i < total_points; i++ )
  {
    C_all.block<3,1>(0,i)  = adapter.getCamOffset(indices[i]);
    X_all.block<3,1>(0,i)  = adapter.getPoint(indices[i]);
    V_all.block<3,1>(0,i)  = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
    A.block<3,1>(3 * i, i) = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
    A.block<3,3>(3 * i, n) = -Eigen::Matrix3d::Identity(3,3);
  }
  //std::cout << "A matrix: " << std::endl << A << std::endl;
  Eigen::MatrixXd U = ((A.transpose() * A).inverse() * A.transpose() ).block(0,0, n, 3 * n);
  //std::cout << "U" << std::endl << U << std::endl;
  /*std::cout << "Data:" << std::endl;
  std::cout << "ci" << std::endl << C_all << std::endl;
  std::cout << "xi" << std::endl << X_all << std::endl;
  std::cout << "vi" << std::endl << V_all << std::endl;*/

  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3,9);
  Eigen::MatrixXd v = Eigen::MatrixXd::Zero(3,1);

  Eigen::MatrixXd fi = Eigen::MatrixXd::Zero(3,1);
  Eigen::MatrixXd pi = Eigen::MatrixXd::Zero(3,1);
  Eigen::MatrixXd vi = Eigen::MatrixXd::Zero(3,1);

  Eigen::MatrixXd fj      = Eigen::MatrixXd::Zero(3,1);
  Eigen::MatrixXd pj      = Eigen::MatrixXd::Zero(3,1);
  Eigen::MatrixXd vj      = Eigen::MatrixXd::Zero(3,1);
  Eigen::MatrixXd u_total = Eigen::MatrixXd::Zero(1, n);
  Eigen::MatrixXd uij     = Eigen::MatrixXd::Zero(3,1);
  
  //Start building the D matrix and v vector from wich the residual can be calculated
  for(int i = 0; i < n; ++i) {
    
    fi = V_all.block<3,1>(0, i);
    pi = X_all.block<3,1>(0, i);
    vi = C_all.block<3,1>(0, i);
    u_total = U.block(i,0, 1, 3 * n);
    double internal_coefficient = 0;
    Eigen::MatrixXd internal_vector = Eigen::MatrixXd::Zero(1,9);
    D = Eigen::MatrixXd::Zero(3,9);
    v = Eigen::MatrixXd::Zero(3,1);
    for(int j = 0; j < n; ++j){
      uij = u_total.block(0, 3 * j , 1, 3);
      pj = X_all.block<3,1>(0, j);
      vj = C_all.block<3,1>(0, j);
      internal_vector.block<1,3>(0,0) = pj(0,0) * uij;
      internal_vector.block<1,3>(0,3) = pj(1,0) * uij;
      internal_vector.block<1,3>(0,6) = pj(2,0) * uij;
     
      internal_coefficient = internal_coefficient + (uij * vj)(0,0);
     
      D = D + (fi * internal_vector);
    }
    
    Eigen::MatrixXd Dr_i = Eigen::MatrixXd::Zero(3,9);
    Dr_i.block<3,3>(0,0) = pi(0,0) * Eigen::MatrixXd::Identity(3,3);
    Dr_i.block<3,3>(0,3) = pi(1,0) * Eigen::MatrixXd::Identity(3,3);
    Dr_i.block<3,3>(0,6) = pi(2,0) * Eigen::MatrixXd::Identity(3,3);
    D = D - Dr_i;
    v = vi - fi * internal_coefficient;
    
    Mr  = Mr  + D.transpose() * D;
    Mrt = Mrt - 2 * D.transpose();
    vr  = vr  + 2 * D.transpose() * v;
    vt  = vt  - 2 * v;
    constant = constant + v.transpose() * v;
  }
  
  Eigen::MatrixXd r = Eigen::MatrixXd::Zero(9,1);
  r(0,0) = rot(0,0);
  r(1,0) = rot(1,0);
  r(2,0) = rot(2,0);
  r(3,0) = rot(0,1);
  r(4,0) = rot(1,1);
  r(5,0) = rot(2,1);
  r(6,0) = rot(0,2);
  r(7,0) = rot(1,2);
  r(8,0) = rot(2,2);
 
  //std::cout << "The squared residual: " << std::endl;
  //std::cout << (r.transpose() * Mr * r + vr.transpose() * r + r.transpose() * Mrt * trans + vt.transpose() * trans + constant + n * trans.transpose() * trans) << std::endl;
  /*std::cout << "Mr:  "      << std::endl << Mr       << std::endl;
  std::cout << "vr:  "      << std::endl << vr       << std::endl;
  std::cout << "Mrt: "      << std::endl << Mrt      << std::endl;
  std::cout << "vt:  "      << std::endl << vt       << std::endl;
  std::cout << "n:   "      << std::endl << n        << std::endl;
  std::cout << "constant: " << std::endl << constant << std::endl;
  std::cout << "R_: "       << std::endl << rot      << std::endl;
  std::cout << "t_: "       << std::endl << trans    << std::endl;*/ 
}

OptimalUPnPFunctionInfo::~OptimalUPnPFunctionInfo(){};

double OptimalUPnPFunctionInfo::objective_function_value(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  const double * p = &rotation(0);
  Map<const Matrix<double,1,9> > r(p, 1, 9);
  Eigen::MatrixXd e = (r * Mr * r.transpose() + vr.transpose() * r.transpose() + r * Mrt * translation + vt.transpose() * translation + n * translation.transpose() * translation);
  return ( e(0,0) );
}

opengv::rotation_t OptimalUPnPFunctionInfo::rotation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  const double * p = &rotation(0);
  Map<const Matrix<double,1,9> > r(p, 1, 9);
  Eigen::MatrixXd result = (2 * Mr * r.transpose()) + ( Mrt * translation ) + vr;
  double * ptr = &result(0);
  Map<Matrix<double, 3,3> > m(ptr, 3, 3);
  return m ;
}

opengv::translation_t OptimalUPnPFunctionInfo::translation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  const double * p = &rotation(0);
  Map<const Matrix<double,1,9> > r(p, 1, 9);
  return (  2 * n * translation + Mrt.transpose() * r.transpose() + vt )  ;
}
