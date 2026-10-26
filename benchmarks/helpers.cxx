
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
#endif

void
timing_start (sc_flopinfo_t &fi, sc_flopinfo_t &snapshot)
{
  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
}

void
timing_end (const char* name, sc_flopinfo_t &fi, sc_flopinfo_t &snapshot, sc_statinfo_t *stats)
{
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, name);
  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
}

void
timing_print_all (sc_statinfo_t *stats, int nstats)
{
  for (int i = 0; i < nstats; i++) {
    sc_stats_print (-1, SC_LP_INFO, 1, &(stats[i]), 1, 1);
  }
}

/**
 * @brief Rank-0-only stdout
 *
 * Every printed value in the examples is global (global element counts,
 * parameters), so the other ranks would only duplicate the lines.
 */
struct root_ostream
{
  bool root;

  template <typename T>
  const root_ostream &
  operator<< (const T &v) const
  {
    if (root)
      std::cout << v;
    return *this;
  }
};

static root_ostream
root_out ()
{
  int rank;
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &rank);
  return { rank == 0 };
}

//=============================================================================
// Test Functions
//=============================================================================

/**
 * @brief Gaussian bump function
 */
template <int U>
auto
gaussian_bump ()
{
  return [] (double x, double y) -> std::array<double, U> {
    const double cx = 0.5, cy = 0.5;
    const double sigma = 0.15;
    const double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
    return { std::exp (-r2 / (2 * sigma * sigma)) };
  };
}

/**
 * @brief Sine wave function
 */
template <int U>
auto
sine_wave ()
{
  return
    [] (double x, double y) -> std::array<double, U> { return { std::sin (4 * M_PI * x) * std::sin (4 * M_PI * y) }; };
}

/**
 * @brief Step function (discontinuous)
 */
template <int U>
auto
step_function ()
{
  return [] (double x, double y) -> std::array<double, U> { return { (x > 0.5 && y > 0.5) ? 1.0 : 0.0 }; };
}

/**
 * @brief Two polynomials separated by a jump along a quarter circle
 */
template <int U>
auto
quarter_circle ()
{
  return [] (double x, double y) -> std::array<double, U> {
    double r = x * x + y * y;
    return { (r < 0.25) ? (x * y + x + 3.) : (x * x * y - 2. * x * y * y + 3. * x) };
  };
}

/**
 * @brief Two components with quarter-circle jumps around the bottom-left
 * (u0) and top-right (u1) corner
 */
template <int U>
auto
two_quarter_circles ()
{
  return [] (double x, double y) -> std::array<double, U> {
    const double xm = 1. - x;
    const double ym = 1. - y;
    const double r0 = x * x + y * y;
    const double r1 = xm * xm + ym * ym;
    return { (r0 < 0.25) ? (x * y + x + 3.) : (x * x * y - 2. * x * y * y + 3. * x),
             (r1 < 0.25) ? (xm * ym + xm + 3.) : (xm * xm * ym - 2. * xm * ym * ym + 3. * xm) };
  };
}

/**
 * @brief 3D Gaussian bump function
 */
template <int U>
auto
gaussian_bump_3d ()
{
  return [] (double x, double y, double z) -> std::array<double, U> {
    const double cx = 0.5, cy = 0.5, cz = 0.5;
    const double sigma = 0.15;
    const double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz);
    return { std::exp (-r2 / (2 * sigma * sigma)) };
  };
}

/**
 * @brief 3D Gaussian bump function
 */
template <int U>
auto
gaussian_bump_3d_periodic ()
{
  return [] (double x, double y, double z) -> std::array<double, U> {
    const double cx = 0.5, cy = 0.5, cz = 0.5;
    const double sigma = 0.15;

    const double xx = x - std::floor(x);
    const double yy = y - std::floor(y);
    const double zz = z - std::floor(z);

    const double r2 = (xx - cx) * (xx - cx) + (yy - cy) * (yy - cy) + (zz - cz) * (zz - cz);
    return { std::exp (-r2 / (2 * sigma * sigma)) };
  };
}

//=============================================================================
// Helpers
//=============================================================================

/**
 * @brief Write high-order Lagrange VTK output (polynomial degree P-1)
 */
template <typename MRA>
void
write_vtk_output (MRA &mra, const std::string &filename)
{

  if (not vtk_activated)
    return;

  root_out () << "  Writing VTK: " << filename << ".vtu\n";

  t8_mra::write_forest_lagrange_vtk (mra, filename.c_str (), MRA::P_DIM - 1);
}

/**
 * @brief Print element and DOF count of the current grid
 */
template <typename MRA>
t8_gloidx_t
print_grid_stats (MRA &mra, const std::string &label)
{
  const auto num_elements = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  const auto num_local_elements = t8_forest_get_local_num_leaf_elements (mra.get_forest ());
  root_out () << "  " << label << ": " << num_local_elements << " local elements, "  << num_elements << " global elements, " << (num_elements * MRA::DOF * MRA::U_DIM)
              << " DOF\n";
  // const auto num_local_elements = t8_forest_get_local_num_leaf_elements (mra.get_forest ());
  // int rank;
  // sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &rank);
  // std::cout << "Process " << rank << " has " << num_local_elements << " local elements" << std::endl;
  return num_elements;
}


/**
 * @brief 3D Gaussian bump function
 */
double
gaussian_bump_3d_periodic_no_mra (double x, double y, double z)
{
    constexpr double cx = 0.5;
    constexpr double cy = 0.5;
    constexpr double cz = 0.5;
    constexpr double sigma = 0.15;

    const double dx = x - cx - std::round(x - cx);
    const double dy = y - cy - std::round(y - cy);
    const double dz = z - cz - std::round(z - cz);

    const double r2 = dx*dx + dy*dy + dz*dz;

    return std::exp(-r2 / (2.0 * sigma * sigma));
}
