#ifndef _ASSEMBLY_H__
#define _ASSEMBLY_H__

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h> // for block structuring
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h> // for block structuring
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h> // Lagrange dg fe elements
//#include <deal.II/fe/fe_dgp.h> // Legendre dg fe elements
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>


namespace Assembly
{



using namespace dealii;

/////////////////////////////////////////////////////////////////////////////
/// Scratch objects for assembling cell matrices and vectors in parallel
/////////////////////////////////////////////////////////////////////////////

/** This is the scratch object which can be used to assemble the local
* cell matrices and vectors.  We do this so that it can be called by
* Workstreams which will distribute the work amoungst all the cores
* by thread building blocks.
*
*/
template<int dim>
struct AssemblyScratch
{
    AssemblyScratch(const FiniteElement<dim> & fe,
                    const Quadrature<dim>		& quadrature,
                    const Quadrature<dim-1>	& face_quadrature);

    AssemblyScratch(const AssemblyScratch & scratch);

    FEValues<dim>					fe_values;
    FEFaceValues<dim>			fe_face_values;

    std::vector<double>					rhs_values;
    std::vector<double>					bc_values;

};

/////////////////////////////////////////////////////////////////////////////
// COPY DATA
////////////////////////////////////////////////////////////////////////////

// Poisson copy data
template<int dim>
struct CopyData
{
    CopyData(const FiniteElement<dim> & fe);

    CopyData(const CopyData & data);

    Vector<double>												local_rhs;
    FullMatrix<double>										local_matrix;
    std::vector<types::global_dof_index>	local_dof_indices;

};


} // end namespace Assembly
#endif
