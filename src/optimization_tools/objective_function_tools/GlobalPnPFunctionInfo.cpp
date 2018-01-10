#include <opengv/optimization_tools/objective_function_tools/GlobalPnPFunctionInfo.hpp>
#include <Eigen/Dense>
#include <opengv/Indices.hpp>
#include <iostream>

GlobalPnPFunctionInfo::GlobalPnPFunctionInfo(const opengv::absolute_pose::AbsoluteAdapterBase & adapter, const opengv::rotation_t & rot, const opengv::translation_t & trans ){

  opengv::Indices indices(adapter.getNumberCorrespondences());
  int total_points = (int) indices.size();
  //std::cout << "Total points: " << total_points << std::endl;
  Mt =  Eigen::Matrix3d::Zero(3,3);
  Mrt = Eigen::MatrixXd::Zero(3,9);
  Mr =  Eigen::MatrixXd::Zero(9,9);
  vt =  Eigen::VectorXd::Zero(3,1);
  vr =  Eigen::VectorXd::Zero(9,1);
  //Eigen::MatrixXd Constant = Eigen::MatrixXd::Zero(1,1);
  
  //Variables defined for debugging
  Eigen::MatrixXd C_all = Eigen::MatrixXd::Zero(3, total_points);
  Eigen::MatrixXd V_all = Eigen::MatrixXd::Zero(3, total_points);
  Eigen::MatrixXd X_all = Eigen::MatrixXd::Zero(3, total_points);
  Eigen::Matrix3d id = Eigen::Matrix3d::Identity(3,3);

  Eigen::Matrix<double,3,1> vi;
  Eigen::Matrix<double,3,1> xi;
  Eigen::Matrix<double,3,1> ci;
  Eigen::Matrix3d Qi;
  Eigen::Matrix3d Vi;
  Eigen::MatrixXd Mr_i  = Eigen::MatrixXd::Zero(3,9);
  Eigen::MatrixXd Mrt_i = Eigen::MatrixXd::Zero(3,9);
  Eigen::Matrix3d I_V   = Eigen::MatrixXd::Zero(3,3);
  Eigen::MatrixXd vr_1  = Eigen::MatrixXd::Zero(1,3);
  Eigen::MatrixXd vr_2  = Eigen::MatrixXd::Zero(1,9);
  for( int i = 0; i < total_points; i++ )
  {
    vi = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
    xi = adapter.getPoint(indices[i]);
    ci = adapter.getCamOffset(indices[i]);
    //C_all.block<3,1>(0,i) = ci;
    //X_all.block<3,1>(0,i) = xi;
    //V_all.block<3,1>(0,i) = vi;
    Vi = vi * vi.transpose() / (vi.transpose() * vi);
    Qi = (id - Vi).transpose() * (id - Vi);
    
    I_V = id - Vi;
    Mt = Mt + Qi;
    vt = vt - (2 * Qi * ci);
    //Constant = Constant + ci.transpose() * Qi * ci;
   
    //Calculate Mr  
    Mr_i.block<3,3>(0,0) = xi(0,0) * I_V;
    Mr_i.block<3,3>(0,3) = xi(1,0) * I_V;
    Mr_i.block<3,3>(0,6) = xi(2,0) * I_V;
   
    Mr = Mr + Mr_i.transpose() * Mr_i;
   
    //Calculate Mrt
    Mrt_i.block<3,3>(0,0) = xi(0,0) * Qi;
    Mrt_i.block<3,3>(0,3) = xi(1,0) * Qi;
    Mrt_i.block<3,3>(0,6) = xi(2,0) * Qi;
    
    Mrt = Mrt + 2*Mrt_i;
       
    //Calculate vr
    vr_1 = ci.transpose() * Qi;
    vr_2.block<1,3>(0,0) = -2 * xi(0,0) * vr_1;
    vr_2.block<1,3>(0,3) = -2 * xi(1,0) * vr_1;
    vr_2.block<1,3>(0,6) = -2 * xi(2,0) * vr_1;
    vr = vr + vr_2.transpose();
  }
  //std::cout << "Data:" << std::endl;
  //std::cout << "ci" << std::endl << C_all << std::endl;
  //std::cout << "xi" << std::endl << X_all << std::endl;
  //std::cout << "vi" << std::endl << V_all << std::endl;
  
  //std::cout << "Final results: " << std::endl;
  //std::cout << "Mt: "  << std::endl << Mt << std::endl;
  //std::cout << "Mrt: " << std::endl << Mrt << std::endl;
  //std::cout << "Mr:"   << std::endl << Mr  << std::endl;
  //std::cout << "vt:"   << std::endl << vt  << std::endl;
  //std::cout << "vr:"   << std::endl << vr   << std::endl;
  //std::cout << "Const: " << std::endl << Constant << std::endl;
}

GlobalPnPFunctionInfo::~GlobalPnPFunctionInfo(){};

double GlobalPnPFunctionInfo::objective_function_value(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  const double * p = &rotation(0);
  Map<const Matrix<double,1,9> > r(p, 1, 9);
  Eigen::MatrixXd e = (translation.transpose() * Mt * translation)
    + (translation.transpose() * Mrt *  r.transpose())
    + (vt.transpose() * translation)
    + (r * Mr * r.transpose() )
    + (vr.transpose() * r.transpose());
  return ( e(0,0));
}

opengv::rotation_t GlobalPnPFunctionInfo::rotation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  const double * p = &rotation(0);
  Map<const Matrix<double,1,9> > r(p, 1, 9);
  Eigen::MatrixXd result = (2 * Mr * r.transpose()) + ( Mrt.transpose() * translation ) + vr;
  double * ptr = &result(0);
  Map<Matrix<double, 3,3> > m(ptr, 3, 3);
  return m;
}

opengv::translation_t GlobalPnPFunctionInfo::translation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  const double * p = &rotation(0);
  Map<const Matrix<double,1,9> > r(p, 1, 9);
  return ( (2 * Mt * translation) + (  Mrt * r.transpose() ) + vt );
}
