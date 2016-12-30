#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h> // for block structuring
#include <deal.II/lac/block_vector.h> // for block structuring
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h> // dg fe elements
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <functional>

// Matrix object which acts as the inverse of a matrix by
// calling an iterative solver.
#include <deal.II/lac/iterative_inverse.h>

// Raviart-Thomas elements
#include <deal.II/fe/fe_raviart_thomas.h>

// Tensor valued functions
#include <deal.II/base/tensor_function.h>

// Multithreading
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>

#include "../source/Assembly.cpp"
#include "../source/functions.cpp"
#include "../source/Grid.cpp"


namespace MixedFEM
{
using namespace dealii;
using namespace std;

template <int dim>
class MixedPoissonProblem
{
public:
    MixedPoissonProblem(const unsigned int degree,
                        const unsigned int n_global_refine,
                        const unsigned int n_local_refine);
    void run();
    void run_with_errors(double & h,
                         double & pot_errors,
                         double & vec_field_errors);


private:
    void make_dofs();
//			void make_Neumann_boundaries();
    void enforce_Neumann_boundary_values();

    void serial_assemble_matrix();
    void serial_assemble_rhs();

    void parallel_assemble_matrix();
    void parallel_assemble_rhs();

    void output_system() const;
    void solve();
    void compute_errors(double & pot_l2_error,
                        double & vec_field_l2_error);
    void output_results() const;

    const unsigned int degree;
    int global_refinement_number;
		int local_refinement_number;

    Triangulation<dim>	triangulation;
    FESystem<dim>		fe;
    DoFHandler<dim>		dof_handler;

    ConstraintMatrix	constraints;

    BlockSparsityPattern		sparsity_pattern;
    BlockSparseMatrix<double>	system_matrix;

    BlockVector<double>			solution;
    BlockVector<double>			system_rhs;

    void assemble_local_matrix(const typename DoFHandler<dim>::active_cell_iterator & cell,
                               Assembly::AssemblyScratch<dim>								& scratch,
                               Assembly::CopyData<dim>											& data,
                               const double 																& lambda);

    void assemble_local_rhs(const typename DoFHandler<dim>::active_cell_iterator & cell,
                            Assembly::AssemblyScratch<dim>								& scratch,
                            Assembly::CopyData<dim>												& data);

    void copy_local_to_global_matrix(const Assembly::CopyData<dim>		& data);
    void copy_local_to_global_rhs(const Assembly::CopyData<dim>				& data);


    SparseDirectUMFPACK direct_solver;
};

}


