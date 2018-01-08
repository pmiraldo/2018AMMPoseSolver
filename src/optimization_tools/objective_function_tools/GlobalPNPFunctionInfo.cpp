#include <opengv/optimization_tools/GlobalPNPFunctionInfo.hpp>
#include <Eigen/Dense>

GlobalPNPFunctionInfo:: GlobalPNPFunctionInfo(const opengv::absolute_pose::AbsoluteAdapterBase & adapter ){

  Eigen::Matrixd3 Mt = Eigen::Matrix3d::Zero(3,3);
  Eigen::MatrixXd Mrt = Eigen::MatrixXd::Zero(9,3);
  Eigen::MatrixXd Mr = Eigen::MatrixXd::Zero(9,9);

  Eigen::Vector3d vt = Eigen::Vector3d::Zero(3,1);
  Eigen::VectorXd vr = Eigen::VectorXd::Zero(9,1);

  Eigen::Matrix3d id = Eigen::Matrix3d::Identity(3,3);

  Eigen::Matrix<double,3,1> vi;
  Eigen::Matrix<double,3,1> xi;
  Eigen::Matrix<double,3,1> ci;
  Eigen::Matrix3d Qi;
  Eigen::Matrix3d Vi;
  Eigen::MatrixXd Mrt_i = Eigen::MatrixXd::Zero(3,9);
  Eigen::MatrixXd Mr1_i = Eigen::MatrixXd::Zero(3,9);
  Eigen::VectorXd Mr2_i = Eigen::VectorXd::Zero(1,9);
  for( int i = 0; i < (int) indices.size(); i++ )
  {
    vi = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
    xi = adapter.getPoint(indices[i]);
    ci = adapter.getCamOffset(indices[i]);

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

    Mr = Mr + (Mr?)
  }

};
GlobalPNPFunctionInfo:: ~GlobalPNPFunctionInfo(){};
