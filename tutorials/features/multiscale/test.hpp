/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2023 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* This file is part of the t8code tutorial regarding the creation of a scmesh.
 * Here, we only wrap the main function. See t8_tutorial_build_cmesh.cxx for the
 * documentation of the 'build-cmesh'-tutorial
 */
#include <t8.h>
#include <t8_forest/t8_forest_general.h>

#include <iostream>
#include <t8_cmesh.hxx>
#include <t8_forest/t8_adapt/multiscale.hpp>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>

void output_test(const t8_mra::t8_multiscale<T8_ECLASS_TRIANGLE>& mst) {
  std::cout << "poly_degree: " << mst.polynomial_degree << "\n";
  std::cout << "dof: " << mst.dof << "\n";
  std::cout << "max_level: " << mst.max_level << "\n";
  std::cout << "c_thres: " << mst.c_thresh << "\n";

  for (auto i = 0u; i < mst.mask_coeffs.size(); ++i) {
    for (auto j = 0u; j < mst.mask_coeffs[i].cols(); ++j) {
      for (auto k = 0u; k < mst.mask_coeffs[i].rows(); ++k)
        std::cout << mst.mask_coeffs[i](k, j) << " ";
      std::cout << "\n";
    }
    std::cout << "\n\n";
  }

  for (auto i = 0u; i < mst.inv_mask_coeffs.size(); ++i) {
    for (auto j = 0u; j < mst.inv_mask_coeffs[i].cols(); ++j) {
      for (auto k = 0u; k < mst.inv_mask_coeffs[i].rows(); ++k)
        std::cout << mst.inv_mask_coeffs[i](k, j) << " ";
      std::cout << "\n";
    }
    std::cout << "\n\n";
  }
}

/* This cmesh is [0,1]² with 8 triangles. */
t8_cmesh_t t8_cmesh_new_basic(sc_MPI_Comm comm) {
  double vertices[72] = {
      0,   0,   0, 0.5, 0,   0, 0.5, 0.5, 0,  // triangle 2
      0,   0,   0, 0,   0.5, 0, 0.5, 0.5, 0,  // triangle 1
      0.5, 0,   0, 1,   0,   0, 1,   0.5, 0,  // triangle 3
      0.5, 0,   0, 0.5, 0.5, 0, 1,   0.5, 0,  // triangle 4
      0,   0.5, 0, 0.5, 0.5, 0, 0.5, 1,   0,  // triangle 5
      0,   0.5, 0, 0,   1,   0, 0.5, 1,   0,  // triangle 6
      0.5, 0.5, 0, 1,   0.5, 0, 1,   1,   0,  // triangle 7
      0.5, 0.5, 0, 0.5, 1,   0, 1,   1,   0,  // triangle 8
  };

  t8_cmesh_t cmesh;
  t8_cmesh_init(&cmesh);

  t8_cmesh_register_geometry<t8_geometry_linear>(cmesh, 2);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class(cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class(cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class(cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class(cmesh, 3, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class(cmesh, 4, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class(cmesh, 5, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class(cmesh, 6, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class(cmesh, 7, T8_ECLASS_TRIANGLE);

  t8_cmesh_set_tree_vertices(cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices(cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices(cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices(cmesh, 3, vertices + 27, 3);
  t8_cmesh_set_tree_vertices(cmesh, 4, vertices + 36, 3);
  t8_cmesh_set_tree_vertices(cmesh, 5, vertices + 45, 3);
  t8_cmesh_set_tree_vertices(cmesh, 6, vertices + 54, 3);
  t8_cmesh_set_tree_vertices(cmesh, 7, vertices + 63, 3);

  t8_cmesh_set_join(cmesh, 0, 1, 1, 1, 0);
  t8_cmesh_set_join(cmesh, 0, 3, 0, 2, 0);
  t8_cmesh_set_join(cmesh, 1, 4, 0, 2, 0);
  t8_cmesh_set_join(cmesh, 2, 3, 1, 1, 0);
  t8_cmesh_set_join(cmesh, 3, 6, 0, 2, 0);
  t8_cmesh_set_join(cmesh, 4, 5, 1, 1, 0);
  t8_cmesh_set_join(cmesh, 4, 7, 0, 2, 0);
  t8_cmesh_set_join(cmesh, 6, 7, 1, 1, 0);

  t8_cmesh_commit(cmesh, comm);

  return cmesh;
}

/// Only for testing
int t8_tutorial_multiscale_main(int argc, char** argv) {
  SC_CHECK_MPI(sc_MPI_Init(&argc, &argv));

  sc_init(sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);
  t8_init(SC_LP_PRODUCTION);
  // sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  t8_mra::t8_multiscale<T8_ECLASS_TRIANGLE> multiscale(3, 0.1, 6u);
  output_test(multiscale);

  // t8_cmesh_t cmesh = example_mesh(comm);
  // t8_scheme_cxx_t* scheme = t8_scheme_new_default_cxx();

  sc_finalize();

  SC_CHECK_MPI(sc_MPI_Finalize());

  return 0;
}
