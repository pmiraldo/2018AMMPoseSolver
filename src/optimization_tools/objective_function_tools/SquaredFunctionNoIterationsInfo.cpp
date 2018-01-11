#include <opengv/optimization_tools/objective_function_tools/SquaredFunctionNoIterationsInfo.hpp>
#include <opengv/Indices.hpp>
#include <iostream>

SquaredFunctionNoIterationsInfo::SquaredFunctionNoIterationsInfo(const opengv::relative_pose::RelativeAdapterBase & adapter ){

  opengv::Indices idx(adapter.getNumberCorrespondences());
   size_t numberCorrespondences = idx.size();
   Eigen::MatrixXd vector_data = Eigen::MatrixXd::Zero (18, 1);
   M = Eigen::MatrixXd::Zero(18,18);
   for( size_t i = 0; i < numberCorrespondences; i++ )
   {
     //std::cout << "iteration: " << i << std::endl;
     opengv::bearingVector_t d1 = adapter.getBearingVector1(idx[i]);
     opengv::bearingVector_t d2 = adapter.getBearingVector2(idx[i]);
     opengv::translation_t v1 = adapter.getCamOffset1(idx[i]);
     opengv::translation_t v2 = adapter.getCamOffset2(idx[i]);
     opengv::rotation_t R1 = adapter.getCamRotation1(idx[i]);
     opengv::rotation_t R2 = adapter.getCamRotation2(idx[i]);

     //unrotate the bearing-vectors to express everything in the body frame
     d1 = R1*d1;
     d2 = R2*d2;

     //generate the Plücker line coordinates
     Eigen::Matrix<double,6,1> l1;
     l1.block<3,1>(0,0) = d1;
     l1.block<3,1>(3,0) = v1.cross(d1);
     Eigen::Matrix<double,6,1> l2;
     l2.block<3,1>(0,0) = d2;
     l2.block<3,1>(3,0) = v2.cross(d2);
     //Calculate the kronecker product
     double a00 = l2(0,0) * l1(0,0); double a01 = l2(0,0) * l1(1,0); double a02 = l2(0,0) * l1(2,0);//e1->First column of the essential matrix
     double a03 = l2(0,0) * l1(3,0); double a04 = l2(0,0) * l1(4,0); double a05 = l2(0,0) * l1(5,0);//r1
     double a06 = l2(1,0) * l1(0,0); double a07 = l2(1,0) * l1(1,0); double a08 = l2(1,0) * l1(2,0);//e2
     double a09 = l2(1,0) * l1(3,0); double a10 = l2(1,0) * l1(4,0); double a11 = l2(1,0) * l1(5,0);//r2
     double a12 = l2(2,0) * l1(0,0); double a13 = l2(2,0) * l1(1,0); double a14 = l2(2,0) * l1(2,0);//e3
     double a15 = l2(2,0) * l1(3,0); double a16 = l2(2,0) * l1(4,0); double a17 = l2(2,0) * l1(5,0);//r3
     double a18 = l2(3,0) * l1(0,0); double a19 = l2(3,0) * l1(1,0); double a20 = l2(3,0) * l1(2,0); //r1->First column of the rotation matrix
     //double a21 = l2(3,0) * l1(3,0); double a22 = l2(3,0) * l1(4,0); double a23 = l2(3,0) * l1(5,0); //0
     double a24 = l2(4,0) * l1(0,0); double a25 = l2(4,0) * l1(1,0); double a26 = l2(4,0) * l1(2,0);//r2
     //double a27 = l2(4,0) * l1(3,0); double a28 = l2(4,0) * l1(4,0); double a29 = l2(4,0) * l1(5,0);//0
     double a30 = l2(5,0) * l1(0,0); double a31 = l2(5,0) * l1(1,0); double a32 = l2(5,0) * l1(2,0);//r3
     //double a33 = l2(5,0) * l1(3,0); double a34 = l2(5,0) * l1(4,0); double a35 = l2(5,0) * l1(5,0);//0
     vector_data(0,0) = a00; vector_data(1,0) = a01; vector_data(2,0) = a02; //e1->Column
     vector_data(3,0) = a06; vector_data(4,0) = a07; vector_data(5,0) = a08; //e2
     vector_data(6,0) = a12; vector_data(7,0) = a13; vector_data(8,0) = a14; //e3


     vector_data(9,0)  = a03 + a18; vector_data(10,0) = a04 + a19; vector_data(11,0) = a05 + a20; //r1->Column R
     vector_data(12,0) = a09 + a24; vector_data(13,0) = a10 + a25; vector_data(14,0) = a11 + a26;
     vector_data(15,0) = a15 + a30; vector_data(16,0) = a16 + a31; vector_data(17,0) = a17 + a32;
    /*std::cout << "l1: "  << " " << l1(0,0) << " " << l1(1,0) << " " << l1(2,0) << " "
              << l1(3,0) << " " << l1(4,0) << " " << l1(5,0) << std::endl;
    std::cout << "l2: "  << " " << l2(0,0) << " " << l2(1,0) << " " << l2(2,0) << " "
    << l2(3,0) << " " << l2(4,0) << " " << l2(5,0) << std::endl;*/
     M = M + (vector_data * vector_data.transpose());
   }
   std::cout << "M matrix: " << std::endl << M << std::endl;

   v = Eigen::MatrixXd::Zero(18,1);
}

SquaredFunctionNoIterationsInfo::~SquaredFunctionNoIterationsInfo(){}

Eigen::MatrixXd SquaredFunctionNoIterationsInfo::get_M(){
  return M;
}

double SquaredFunctionNoIterationsInfo::objective_function_value(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  //Eigen::MatrixXd v = Eigen::MatrixXd::Zero(18,1);

  v(0,0)  = translation(1,0) * rotation(2,0) - translation(2,0) * rotation(1,0);
  v(1,0)  = translation(2,0) * rotation(0,0) - translation(0,0) * rotation(2,0);
  v(2,0)  = translation(0,0) * rotation(1,0) - translation(1,0) * rotation(0,0);
  v(3,0)  = translation(1,0) * rotation(2,1) - translation(2,0) * rotation(1,1);
  v(4,0)  = translation(2,0) * rotation(0,1) - translation(0,0) * rotation(2,1);
  v(5,0)  = translation(0,0) * rotation(1,1) - translation(1,0) * rotation(0,1);
  v(6,0)  = translation(1,0) * rotation(2,2) - translation(2,0) * rotation(1,2);
  v(7,0)  = translation(2,0) * rotation(0,2) - translation(0,0) * rotation(2,2);
  v(8,0)  = translation(0,0) * rotation(1,2) - translation(1,0) * rotation(0,2);
  
  v(9,0)  = rotation(0,0);
  v(10,0) = rotation(1,0);
  v(11,0) = rotation(2,0);
  v(12,0) = rotation(0,1);
  v(13,0) = rotation(1,1);
  v(14,0) = rotation(2,1);
  v(15,0) = rotation(0,2);
  v(16,0) = rotation(1,2);
  v(17,0) = rotation(2,2);


  return ( (v.transpose() * M * v ) (0,0) );
}

opengv::rotation_t SquaredFunctionNoIterationsInfo::rotation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  //Eigen::MatrixXd v = Eigen::MatrixXd::Zero(18,1);

  v(0,0)  = translation(1,0) * rotation(2,0) - translation(2,0) * rotation(1,0);
  v(1,0)  = translation(2,0) * rotation(0,0) - translation(0,0) * rotation(2,0);
  v(2,0)  = translation(0,0) * rotation(1,0) - translation(1,0) * rotation(0,0);
  v(3,0)  = translation(1,0) * rotation(2,1) - translation(2,0) * rotation(1,1);
  v(4,0)  = translation(2,0) * rotation(0,1) - translation(0,0) * rotation(2,1);
  v(5,0)  = translation(0,0) * rotation(1,1) - translation(1,0) * rotation(0,1);
  v(6,0)  = translation(1,0) * rotation(2,2) - translation(2,0) * rotation(1,2);
  v(7,0)  = translation(2,0) * rotation(0,2) - translation(0,0) * rotation(2,2);
  v(8,0)  = translation(0,0) * rotation(1,2) - translation(1,0) * rotation(0,2);
  
  v(9,0)  = rotation(0,0);
  v(10,0) = rotation(1,0);
  v(11,0) = rotation(2,0);
  v(12,0) = rotation(0,1);
  v(13,0) = rotation(1,1);
  v(14,0) = rotation(2,1);
  v(15,0) = rotation(0,2);
  v(16,0) = rotation(1,2);
  v(17,0) = rotation(2,2);
  v = M * v;
  Eigen::MatrixXd skew_symmetric = Eigen::MatrixXd::Zero(3,3);
  skew_symmetric(0,1) =  translation(2,0); skew_symmetric(0,2) = -translation(1,0);
  skew_symmetric(1,0) = -translation(2,0); skew_symmetric(1,2) =  translation(0,0);
  skew_symmetric(2,0) =  translation(1,0); skew_symmetric(2,1) = -translation(0,0);
  
  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(3,3);
  grad.block<3,1>(0,0) = skew_symmetric * v.block<3,1>(0,0) + v.block<3,1>(9,0);
  grad.block<3,1>(0,1) = skew_symmetric * v.block<3,1>(3,0) + v.block<3,1>(12,0);
  grad.block<3,1>(0,2) = skew_symmetric * v.block<3,1>(6,0) + v.block<3,1>(15,0);
  return 2 * grad;
}

opengv::translation_t SquaredFunctionNoIterationsInfo::translation_gradient(const opengv::rotation_t & rotation, const opengv::translation_t & translation){
  Eigen::MatrixXd dv = Eigen::MatrixXd::Zero(3,9);
 
  dv(0,0) = 0;             dv(0,1) = -rotation(2,0); dv(0,2) = rotation(1,0);
  dv(0,3) = 0;             dv(0,4) = -rotation(2,1); dv(0,5) = rotation(1,1);
  dv(0,6) = 0;             dv(0,7) = -rotation(2,2); dv(0,8) = rotation(1,2);

  dv(1,0) = rotation(2,0); dv(1,1) = 0;              dv(1,2) = -rotation(0,0);
  dv(1,3) = rotation(2,1); dv(1,4) = 0;              dv(1,5) = -rotation(0,1);
  dv(1,6) = rotation(2,2); dv(1,7) = 0;              dv(1,8) = -rotation(0,2);

  dv(2,0) = -rotation(1,0); dv(2,1) = rotation(0,0); dv(2,2) = 0;
  dv(2,3) = -rotation(1,1); dv(2,4) = rotation(0,1); dv(2,5) = 0;
  dv(2,6) = -rotation(1,2); dv(2,7) = rotation(0,2); dv(2,8) = 0;

  Eigen::MatrixXd v = Eigen::MatrixXd::Zero(18,1);

  v(0,0)  = translation(1,0) * rotation(2,0) - translation(2,0) * rotation(1,0);
  v(1,0)  = translation(2,0) * rotation(0,0) - translation(0,0) * rotation(2,0);
  v(2,0)  = translation(0,0) * rotation(1,0) - translation(1,0) * rotation(0,0);
  v(3,0)  = translation(1,0) * rotation(2,1) - translation(2,0) * rotation(1,1);
  v(4,0)  = translation(2,0) * rotation(0,1) - translation(0,0) * rotation(2,1);
  v(5,0)  = translation(0,0) * rotation(1,1) - translation(1,0) * rotation(0,1);
  v(6,0)  = translation(1,0) * rotation(2,2) - translation(2,0) * rotation(1,2);
  v(7,0)  = translation(2,0) * rotation(0,2) - translation(0,0) * rotation(2,2);
  v(8,0)  = translation(0,0) * rotation(1,2) - translation(1,0) * rotation(0,2);
  
  v(9,0)  = rotation(0,0);
  v(10,0) = rotation(1,0);
  v(11,0) = rotation(2,0);
  v(12,0) = rotation(0,1);
  v(13,0) = rotation(1,1);
  v(14,0) = rotation(2,1);
  v(15,0) = rotation(0,2);
  v(16,0) = rotation(1,2);
  v(17,0) = rotation(2,2);


  return (2 * dv * M.block<9,18>(0,0) * v);
}
