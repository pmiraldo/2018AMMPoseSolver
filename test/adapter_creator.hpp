#ifndef ADAPTER_CREATOR_HPP
#define ADAPTER_CREATOR_HPP
#include <opengv/absolute_pose/NoncentralAbsoluteAdapter.hpp>
#include <opengv/types.hpp>

void create_random_pose(opengv::rotation_t & rotation, opengv::translation_t & position, const double & noise, const double & outlierFraction, const size_t & numberPoints, const int & numberCameras, opengv::bearingVectors_t & bearingVectors, std::vector<int> & camCorrespondences, opengv::points_t & points, opengv::translations_t & camOffsets, opengv::rotations_t & camRotations);



void choose_best_transformation(const opengv::transformations_t &solutions, const opengv::transformation_t & ref, opengv::transformation_t & transformation_method);

void choose_best_essential_matrix(const opengv::essentials_t & solutions, const opengv::essential_t & ref, opengv::essential_t & essential_matrix);

opengv::essential_t calculate_essential_matrix(const opengv::translation_t & position, const opengv::rotation_t & rotation);

opengv::transformation_t calculate_transformation_from_essential_matrix(const opengv::essential_t & essential_matrix, const opengv::rotation_t &ref_rotation, const opengv::translation_t & ref_translation);

#endif
