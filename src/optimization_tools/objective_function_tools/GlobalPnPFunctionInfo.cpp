#include <opengv/optimization_tools/objective_function_tools/GlobalPnPFunctionInfo.hpp>
#include <Eigen/Dense>
#include <opengv/Indices.hpp>
#include <iostream>

GlobalPNPFunctionInfo:: GlobalPNPFunctionInfo(const opengv::absolute_pose::AbsoluteAdapterBase & adapter ){

  opengv::Indices indices(adapter.getNumberCorrespondences());
  int total_points = (int) indices.size();
  Eigen::Matrix3d Mt = Eigen::Matrix3d::Zero(3,3);
  Eigen::MatrixXd Mrt = Eigen::MatrixXd::Zero(9,3);
  Eigen::MatrixXd Mr = Eigen::MatrixXd::Zero(9,9);
  Eigen::Vector3d vt = Eigen::Vector3d::Zero(3,1);
  Eigen::VectorXd vr = Eigen::VectorXd::Zero(9,1);
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
  Eigen::MatrixXd Mrt_i = Eigen::MatrixXd::Zero(3,9);
  Eigen::MatrixXd Mr1_i = Eigen::MatrixXd::Zero(3,9);
  Eigen::VectorXd Mr2_i = Eigen::VectorXd::Zero(1,9);
  Eigen::Vector3d vr_1 = Eigen::Vector3d::Zero(1,3);
  for( int i = 0; i < total_points; i++ )
  {
    vi = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
    xi = adapter.getPoint(indices[i]);
    ci = adapter.getCamOffset(indices[i]);
    C_all.block<3,1>(0,i) = ci;
    X_all.block<3,1>(0,i) = xi;
    V_all.block<3,1>(0,i) = vi;
    Vi = vi * vi.transpose() / (vi.transpose() * vi);
    Qi = (id - Vi).transpose() * (id - Vi);
    
    Mt = Mt + Qi;
    vt = vt - (2 * Qi * ci);

    //Calculate Mrt
    Mrt_i.block<1,3>(0,0) = Qi(0,0) * xi.transpose(); Mrt_i.block<1,3>(0,3) = Qi(0,1) * xi.transpose(); Mrt_i.block<1,3>(0,6) = Qi(0,2) * xi.transpose();
    Mrt_i.block<1,3>(1,0) = Qi(1,0) * xi.transpose(); Mrt_i.block<1,3>(1,3) = Qi(1,1) * xi.transpose(); Mrt_i.block<1,3>(1,6) = Qi(1,2) * xi.transpose();
    Mrt_i.block<1,3>(2,0) = Qi(2,0) * xi.transpose(); Mrt_i.block<1,3>(2,3) = Qi(2,1) * xi.transpose(); Mrt_i.block<1,3>(2,6) = Qi(2,2) * xi.transpose();

    Mrt = Mrt + Mrt_i.transpose();

    //Calculate Mr
    Mr1_i.block<1,3>(0,0) = Vi(0,0) * xi.transpose();Mr1_i.block<1,3>(0,3) = Vi(0,1) * xi.transpose();Mr1_i.block<1,3>(0,6) = Vi(0,2) * xi.transpose();
    Mr1_i.block<1,3>(1,0) = Vi(1,0) * xi.transpose();Mr1_i.block<1,3>(1,3) = Vi(1,1) * xi.transpose();Mr1_i.block<1,3>(1,6) = Vi(1,2) * xi.transpose();
    Mr1_i.block<1,3>(2,0) = Vi(2,0) * xi.transpose();Mr1_i.block<1,3>(2,3) = Vi(2,1) * xi.transpose();Mr1_i.block<1,3>(2,6) = Vi(2,2) * xi.transpose();

    Mr2_i.block<1,3>(0,0) = vi(0,0) * xi.transpose();Mr2_i.block<1,3>(1,0) = vi(0,0) * xi.transpose();Mr2_i.block<1,3>(0,0) = vi(2,0) * xi.transpose();
    double alpha  = vi.transpose() * vi;

    Mr = Mr + (  (Mr1_i.transpose() * Mr1_i)  - (2/alpha) * Mr2_i.transpose() * Mr2_i );

    //Calculate vr
    vr_1 = ci.transpose() * Qi;
    vr(0,0) = vr_1(0,0) * xi(0,0);vr(1,0) = vr_1(0,0) * xi(1,0);vr(2,0) = vr_1(0,0) * xi(2,0);
    vr(4,0) = vr_1(0,1) * xi(0,0);vr(5,0) = vr_1(0,1) * xi(1,0);vr(6,0) = vr_1(0,1) * xi(2,0);
    vr(7,0) = vr_1(0,2) * xi(0,0);vr(8,0) = vr_1(0,2) * xi(1,0);vr(9,0) = vr_1(0,2) * xi(2,0);
    
  }
  std::cout << "Data:" << std::endl;
  std::cout << "ci" << std::endl << C_all << std::endl;
  std::cout << "xi" << std::endl << X_all << std::endl;
  std::cout << "vi" << std::endl << V_all << std::endl;

  std::cout << "Final results: " << std::endl;
  std::cout << "Mt: "  << std::endl << Mt << std::endl;
  std::cout << "Mrt: " << std::endl << Mrt << std::endl;
  std::cout << "Mr:"   << std::endl << Mr  << std::endl;
  std::cout << "vt:"   << std::endl << vt  << std::endl;
  std::cout << "vr:"   << std::endl << vr   << std::endl; 
}

GlobalPNPFunctionInfo:: ~GlobalPNPFunctionInfo(){};
