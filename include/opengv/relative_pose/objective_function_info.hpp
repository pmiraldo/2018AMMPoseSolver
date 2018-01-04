#ifndef OBJECTIVE_FUNCTION_INFO_H
#define OBJECTIVE_FUNCTION_INFO_H

#include <Eigen/Dense>
#include <opengv/types.hpp>


class objective_function_info{
 public:

   
  virtual double objective_function_value(const opengv::rotation_t & rotation,
					  const opengv::translation_t & translation) = 0;
  virtual opengv::rotation_t rotation_gradient(const opengv::rotation_t & rotation,
				       const opengv::translation_t & translation) = 0;
  virtual opengv::translation_t translation_gradient(const opengv::rotation_t & rotation,
					     const opengv::translation_t & translation) = 0;
  
  //virtual objective_function_info(Eigen::MatrixXd LPluckerC1, Eigen::MatrixXd LPluckerW2);
  //virtual ~objective_function_info();
  //private:
  //Eigen::MatrixXd M;
};

#endif
