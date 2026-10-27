/**
 * @file t8_time_mra_1
 */
#include "t8.h"
#if T8_ENABLE_MRA

#include "t8_mra/t8_mra.hxx"
#include "t8_cmesh/t8_cmesh.h"
#include "t8_cmesh/t8_cmesh_examples.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

#include <sc_refcount.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

#include <t8_vtk/t8_vtk_writer.h>
#include <t8_types/t8_vec.hxx>

//=============================================================================
// File-wide settings
//=============================================================================

// Handle vtk output
#if T8_ENABLE_DEBUG
  const bool vtk_activated = true;
#else
  const bool vtk_activated = true;
#endif

// To switch between different paramter sets (coarser, finer, etc.)
#define PARAMETER_SET 3

// Lazy include to outsource helper functions
#include "helpers.cxx"

//=============================================================================
// Test 1: Hex (3D), uniform initialization
//=============================================================================

void
example_hex_3d (int n_blocks_per_dir)
{
  root_out () << "\n=== 5. Hex: 3D adaptation ===\n";

  // Some variables for the number of blocks
  const int n_x = n_blocks_per_dir;
  const int n_y = n_x;
  const int n_z = n_x;
  const int num_global_trees = n_x * n_y * n_z;

#if PARAMETER_SET == 1
    // Define parameters used in this example.
    constexpr int U = 1;
    constexpr int P = 3;
    const int min_level = 0;
    const int max_level = 5;
    const double c_thresh = 250.0; /// double(num_global_trees);
    const int gamma = 3;
#elif PARAMETER_SET == 2
    // Define parameters used in this example.
    constexpr int U = 1;
    constexpr int P = 3;
    const int min_level = 0;
    const int max_level = 6;
    const double c_thresh = 300.0; /// double(num_global_trees);
    const int gamma = 3;
#elif PARAMETER_SET == 3
    // Define parameters used in this example.
    constexpr int U = 1;
    constexpr int P = 3;
    const int min_level = 0;
    const int max_level = 7;
    const double c_thresh = 1000.0; /// double(num_global_trees);
    const int gamma = 3;
#endif

  // Set up MRA instance.
  t8_mra::multiscale<T8_ECLASS_HEX, U, P> mra (max_level, sc_MPI_COMM_WORLD);

  // Define mesh and scheme.
  const t8_gloidx_t local_tree_offset = 0;
  const double boundary_coords[24] = { 0, 0, 0, n_x, 0, 0, 0, n_y, 0, n_x, n_y, 0, 0, 0, n_z, n_x, 0, n_z, 0, n_y, n_z, n_x, n_y, n_z };
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_pad_ext (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, boundary_coords, n_x, n_y,
                                                     n_z, 0, 0, 0, false, false, local_tree_offset);
  auto *scheme = t8_scheme_new_default ();

  // Vtk output
  if (vtk_activated)
    t8_cmesh_vtk_write_file (cmesh, "mra_output/mra_cmesh_hexes");

  // Prepare profiling
  sc_flopinfo_t fi, snapshot;
  int const num_stats = 2;
  sc_statinfo_t stats[num_stats];

  // Step 0; Initialize the MRA grid on maximum level.
  timing_start (fi, snapshot);
  mra.initialize_data (cmesh, scheme, max_level, gaussian_bump_3d_periodic<U> ());
  timing_end ("MRA_Init", fi, snapshot, &(stats[0]));

  // Output
  print_grid_stats (mra, "Uniform level " + std::to_string (max_level));
  write_vtk_output (mra, "mra_output/05_hex_step0_uniform");

  // Step 1: Coarsen the grid.
  timing_start (fi, snapshot);
  mra.coarsen (min_level, max_level, t8_mra::hard_thresholding { .c_thresh = c_thresh, .gamma = gamma });
  print_grid_stats (mra, "After coarsening");
  timing_end ("MRA_Coarsen", fi, snapshot, &(stats[1]));
  write_vtk_output (mra, "mra_output/05_hex_step1_coarsened");

  // Print runtime statistics
  timing_print_all(stats,num_stats);

  // Clean up memory.
  mra.cleanup ();

  // Print hint to user.
  root_out () << "  Use ParaView 'Clip' / 'Slice' filters to see the internal structure.\n";
}

//=============================================================================
// Test 2: Hex (3D), adaptive bottom-up initialization
//=============================================================================

void
example_hex_3d_bottom_up (int n_blocks_per_dir)
{
  root_out () << "\n=== 5. Hex: 3D adaptation ===\n";


  // Number of blocks
  const int n_x = n_blocks_per_dir;
  const int n_y = n_x;
  const int n_z = n_x;
  const int num_global_trees = n_x * n_y * n_z;


#if PARAMETER_SET == 1
  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 6;
  const double c_thresh = 0.35; /// double(num_global_trees);
  const int gamma = 1;
#elif PARAMETER_SET == 2
  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 7;
  const double c_thresh = 150.0; /// double(num_global_trees);
  const int gamma = 2;
#elif PARAMETER_SET == 3
  constexpr int U = 1;
  constexpr int P = 3;
  const int min_level = 0;
  const int max_level = 8;
  const double c_thresh = 150.0; /// double(num_global_trees);
  const int gamma = 2;
#endif

  // Prepare profiling
  sc_flopinfo_t fi, snapshot;
  int const num_stats = 4;
  sc_statinfo_t stats[num_stats];

  // Prepare cmesh.
  const t8_gloidx_t local_tree_offset = 0;
  const double boundary_coords[24] = { 0, 0, 0, n_x, 0, 0, 0, n_y, 0, n_x, n_y, 0, 0, 0, n_z, n_x, 0, n_z, 0, n_y, n_z, n_x, n_y, n_z };

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_pad_ext (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, boundary_coords, n_x, n_y,
                                                     n_z, 0, 0, 0, false, false, local_tree_offset);
  // vtk output of cmesh
  if (vtk_activated)
    t8_cmesh_vtk_write_file (cmesh, "mra_output/mra_cmesh_hexes");

  // Iterate over different versions of refinement (refine-only, grade-only, balance)
  for(int i=0; i<3; i++)
  {

    root_out() << "Refine Version Iterator: " << i << "\n";

    /// Increase cmesh ref counter to avoid destruction
    if(i<2) t8_cmesh_ref(cmesh);

    // Set up MRA instance.
    t8_mra::multiscale<T8_ECLASS_HEX, U, P> mra (max_level, sc_MPI_COMM_WORLD);

    // Step 0; Initialize the MRA grid on maximum level.
    timing_start (fi, snapshot);
    mra.initialize_data_adaptive (cmesh, t8_scheme_new_default (), max_level, gaussian_bump_3d_periodic<U> (),
                                  t8_mra::hard_thresholding { .c_thresh = c_thresh, .gamma = gamma});
    timing_end ("MRA_Bottom_Up", fi, snapshot, &(stats[0]));

    // Output
    print_grid_stats (mra, "Bottom Up");
    write_vtk_output (mra, "mra_output/05_hex_step0_bottom_up");

    // Step 1: Refine, different versions.
    if(i==0) {
      timing_start (fi, snapshot);
      mra.refine (min_level, max_level, t8_mra::harten_prediction { .c_thresh = 0.025 * c_thresh, .gamma = gamma, .switch_off_grading = false, .switch_off_refining = true });
      timing_end ("MRA_Grade", fi, snapshot, &(stats[1]));
      print_grid_stats (mra, "After refinement A (only grading)");
      write_vtk_output (mra, "mra_output/05_hex_step1A_graded");
    }
    else if (i==1) {
      timing_start (fi, snapshot);
      mra.refine (min_level, max_level, t8_mra::harten_prediction { .c_thresh = 0.1 * c_thresh, .gamma = gamma, .switch_off_grading = true, .switch_off_refining = false });
      timing_end ("MRA_Refine", fi, snapshot, &(stats[2]));
      print_grid_stats (mra, "After refinement B (no grading)");
      write_vtk_output (mra, "mra_output/05_hex_step1B_refined");
    }
    else {
        timing_start (fi, snapshot);
        mra.balance();
        timing_end ("MRA_Balance", fi, snapshot, &(stats[3]));
        print_grid_stats (mra, "Balanced");
        write_vtk_output (mra, "mra_output/05_hex_step1C_balanced");
    }

    // Clean up memory.
    mra.cleanup ();

  }

  // Print stats
  timing_print_all(stats,num_stats);

  // Print hint to user.
  root_out () << "  Use ParaView 'Clip' / 'Slice' filters to see the internal structure.\n";
}

//=============================================================================
//=============================================================================



int
t8_no_mra_adapt_callback ([[maybe_unused]] t8_forest_t forest,  [[maybe_unused]] t8_forest_t forest_from,
                          [[maybe_unused]] t8_locidx_t which_tree,
                          [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                          [[maybe_unused]] const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                          [[maybe_unused]] const int num_elements, [[maybe_unused]] t8_element_t *elements[])
{

  t8_3D_vec centroid;
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid.data());
  double fun_eval = t8_dist(centroid, t8_3D_vec{0.5, 0.5, 0.5});

  int level = scheme->element_get_level(tree_class, elements[0]);

  if( is_family and level > 1 and fun_eval > 0.5) return -1;
  return 0;
}


//=============================================================================
// Test 3: Hex (3D), no-MRA version
//=============================================================================

void
example_hex_3d_no_mra (int n_blocks_per_dir)
{
  root_out () << "\n=== 5. Hex: 3D adaptation ===\n";

  // Number of blocks.
  const int n_x = n_blocks_per_dir;
  const int n_y = n_x;
  const int n_z = n_x;
  const int num_global_trees = n_x * n_y * n_z;

#if PARAMETER_SET == 1
 const int start_level = 5;
#elif PARAMETER_SET == 2
 const int start_level = 6;
#elif PARAMETER_SET == 3
 const int start_level = 7;
#endif

  // Define mesh and scheme.
  const t8_gloidx_t local_tree_offset = 0;
  const double boundary_coords[24] = { 0, 0, 0, n_x, 0, 0, 0, n_y, 0, n_x, n_y, 0, 0, 0, n_z, n_x, 0, n_z, 0, n_y, n_z, n_x, n_y, n_z };

  // Set cmesh and scheme
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_pad_ext (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, boundary_coords, n_x, n_y,
                                                     n_z, 0, 0, 0, false, false, local_tree_offset);

  auto *scheme = t8_scheme_new_default ();
  if (vtk_activated)
    t8_cmesh_vtk_write_file (cmesh, "no_mra_output/mra_cmesh_hexes");

  // Prepare profiling
  sc_flopinfo_t fi, snapshot;
  int const num_stats = 3;
  sc_statinfo_t stats[num_stats];

  // Step 0; Initialize the MRA grid on maximum level.
  timing_start (fi, snapshot);
  t8_forest_t uniform_forest = t8_forest_new_uniform(cmesh, scheme, start_level, 1, sc_MPI_COMM_WORLD);
  timing_end ("New_Uniform", fi, snapshot, &(stats[0]));
  if (vtk_activated) t8_forest_write_vtk(uniform_forest, "no_mra_output/uniform_forest");

  root_out () << "Uniform forest on level " << start_level << ":    " << t8_forest_get_global_num_leaf_elements(uniform_forest) << " elements\n";

  // Step 1: Adapt the grid.
  timing_start (fi, snapshot);
  t8_forest_t adapted_forest = t8_forest_new_adapt(uniform_forest, t8_no_mra_adapt_callback, 1, 1, nullptr );
  timing_end ("Adapt", fi, snapshot, &(stats[1]));

  root_out () << "Adapted forest:    " << t8_forest_get_global_num_leaf_elements(adapted_forest)  << " elements\n";
  if (vtk_activated) t8_forest_write_vtk(adapted_forest, "no_mra_output/adapted_forest");

  // Step 2: Balance
  t8_forest_t balanced_forest;
  t8_forest_init(&balanced_forest);
  t8_forest_set_balance (balanced_forest, adapted_forest, 0);
  timing_start (fi, snapshot);
  t8_forest_commit(balanced_forest);
  timing_end ("Balance", fi, snapshot, &(stats[2]));

  root_out () << "Balanced forest:    " << t8_forest_get_local_num_leaf_elements(balanced_forest)  << " local elements,    "<< t8_forest_get_global_num_leaf_elements(balanced_forest)  << " global elements\n";

  // Vtk output
  if (vtk_activated) t8_forest_write_vtk(balanced_forest, "no_mra_output/balanced_forest");


  // Print statistics
  timing_print_all(stats,num_stats);

  // Clean up memory
  t8_forest_unref(&balanced_forest);

}


//=============================================================================
// Main
//=============================================================================

/**
 * \brief Main routine running all six examples covered in this file.
 *
 * The examples are:
 *  1. Full adaptation cycle (top-down): initialize -> coarsen -> refine -> coarsen
 *  2. Bottom-up initialization: adaptive grid without building the uniform grid
 *  3. Custom adaptation criteria
 *  4. Triangle vs quad comparison on the same data
 *  5. 3D (hex) adaptation
 *  6. Two state variables (U = 2) with different jump locations
 *
 * \param [in] argc Standard argc (not used, only passed to sc_MPI_Init)
 * \param [in] argv Standard argv (not used, only passed to sc_MPI_Init)
 * \return 0
 */
int
main (int argc, char **argv)
{
  // Initialize SC library for MPI handling.
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, nullptr, SC_LP_ESSENTIAL);

  // Initialize t8code.
  t8_init (SC_LP_ESSENTIAL);

  // Read command-line arguments
  int num_blocks_per_direction = std::atoi(argv[1]);

  // Rank 0 creates the output directory; everyone writes into it
  int mpirank;
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  if (mpirank == 0)
    std::filesystem::create_directory ("mra_output");
  sc_MPI_Barrier (sc_MPI_COMM_WORLD);
  if (!std::filesystem::exists ("mra_output"))
    t8_errorf ("Could not create directory");


  root_out() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  root_out() << "++++++++++++++++ TEST 1: MRA with uniform initialization ++++++++++++++++++++\n";
  root_out() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  example_hex_3d (num_blocks_per_direction);
  root_out() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  root_out() << "++++++++++++++++ TEST 2: MRA with adaptive initialization +++++++++++++++++++\n";
  root_out() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  example_hex_3d_bottom_up (num_blocks_per_direction);
  root_out() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  root_out() << "+++++++++++++++++++++++++++++++ TEST 3: No MRA ++++++++++++++++++++++++++++++\n";
  root_out() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  example_hex_3d_no_mra (num_blocks_per_direction);
  // Log output.
  root_out () << "\nAll examples completed. Output in mra_output/ (open in ParaView).\n";

  // Finalize SC library.
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

#endif  // T8_ENABLE_MRA
