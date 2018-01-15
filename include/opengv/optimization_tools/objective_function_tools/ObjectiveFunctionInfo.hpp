#ifndef OBJECTIVEFUNCTIONINFO_H
#define OBJECTIVEFUNCTIONINFO_H

#include <Eigen/Dense>
#include <opengv/types.hpp>


class ObjectiveFunctionInfo{
 public:

   
  virtual double objective_function_value(const opengv::rotation_t & rotation, const opengv::translation_t & translation) = 0;
  virtual opengv::rotation_t rotation_gradient(const opengv::rotation_t & rotation,
				       const opengv::translation_t & translation) = 0;
  virtual opengv::translation_t translation_gradient(const opengv::rotation_t & rotation,
					     const opengv::translation_t & translation) = 0;
};

#endif
