#include "../include/Assembly.hpp"

namespace Assembly
{

/////////////////////////////////////////////////////////////////////////////
/// Scratch objects for assembling cell matrices and vectors in parallel
/////////////////////////////////////////////////////////////////////////////

// constructor for the Asssembly Scratch
template<int dim>
AssemblyScratch<dim>::
AssemblyScratch(const FiniteElement<dim> & fe,
                const Quadrature<dim>		& quadrature,
                const Quadrature<dim-1>	& face_quadrature)
    :
    fe_values(fe, quadrature,
             update_values 					 |
             update_gradients 			 |
             update_quadrature_points |
             update_JxW_values				 ),
    fe_face_values(fe, face_quadrature,
                  update_values 					 |
                  update_normal_vectors		 |
                  update_quadrature_points |
                  update_JxW_values				 ),
    rhs_values(quadrature.size()),
    bc_values(face_quadrature.size())
{}

/**  The constructor which you will call if you want to create an
* AssemblyScratch object which will be called when you are looping
* over all the cells by hand and assembling the local conributions
* in a sequential manner.
*/


// copy constructor for the Asssembly Scratch
template<int dim>
AssemblyScratch<dim>::
AssemblyScratch(const AssemblyScratch & scratch)
    :
    fe_values(scratch.fe_values.get_fe(),
             scratch.fe_values.get_quadrature(),
             scratch.fe_values.get_update_flags() ),
    fe_face_values(scratch.fe_face_values.get_fe(),
                  scratch.fe_face_values.get_quadrature(),
                  scratch.fe_face_values.get_update_flags() ),
    rhs_values(scratch.rhs_values),
    bc_values(scratch.bc_values)
{}


// constructor
template<int dim>
CopyData<dim>::
CopyData(const FiniteElement<dim> & fe)
    :
    local_rhs(fe.dofs_per_cell),
    local_matrix(fe.dofs_per_cell,
                fe.dofs_per_cell),
    local_dof_indices(fe.dofs_per_cell)
{ }

// copy constructor
template<int dim>
CopyData<dim>::
CopyData(const CopyData & data)
    :
    local_rhs(data.local_rhs),
    local_matrix(data.local_matrix),
    local_dof_indices(data.local_dof_indices)
{ }

}