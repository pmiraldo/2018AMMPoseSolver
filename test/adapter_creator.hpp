#ifndef ADAPTER_CREATOR_HPP
#define ADAPTER_CREATOR_HPP
#include <opengv/absolute_pose/NoncentralAbsoluteAdapter.hpp>
#include <opengv/types.hpp>

void create_random_pose(opengv::rotation_t & rotation, opengv::translation_t & position, const double & noise, const double & outlierFraction, const size_t & numberPoints, const int & numberCameras, opengv::bearingVectors_t & bearingVectors, std::vector<int> & camCorrespondences, opengv::points_t & points, opengv::translations_t & camOffsets, opengv::rotations_t & camRotations);



void choose_best_rotation(const opengv::transformations_t &solutions, const opengv::rotation_t & rotation, opengv::transformation_t & transformation_method);

#endif
