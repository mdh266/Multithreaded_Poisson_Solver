#include "../include/Mixed.hpp"
namespace MixedFEM
{
using namespace dealii;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// MIXED POISSON CLASS
///////////////////////////////////////////////////////////////////////////////
template <int dim>
MixedPoissonProblem<dim>::MixedPoissonProblem(const unsigned int degree,
        const unsigned int n_global_refine,
        const unsigned int n_local_refine)
    :
    degree(degree),
    global_refinement_number(n_global_refine),
    local_refinement_number(n_local_refine),
    fe(FE_RaviartThomas<dim>(degree), 1,
      FE_DGQ<dim>(degree), 1),
    dof_handler(triangulation)
{}


template <int dim>
void MixedPoissonProblem<dim>::make_dofs()
{

    dof_handler.distribute_dofs(fe);

    // Renumber dofs for [ VectorField, Potential ]^{T} set up
    DoFRenumbering::component_wise(dof_handler);

    std::vector<types::global_dof_index> dofs_per_component(dim+1);
    DoFTools::count_dofs_per_component(dof_handler, dofs_per_component);

    // get number of dofs in vector field components and of potential
    // in each component/dimension of vector field has same number of dofs
    // NOTE: Raviart-Thomas elements are not primitive, that is why we skip
    // from component 0 to compoenent dim.
    const unsigned int n_vector_field = dofs_per_component[0];
    const unsigned int n_potential = dofs_per_component[dim];

    std::cout << "Number of active cells : "
              << triangulation.n_active_cells()
              << std::endl
              << "Total number of cells: "
              << triangulation.n_cells()
              << std::endl
              << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << " (" << n_vector_field << " + " << n_potential << ")"
              << std::endl;

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            constraints);

    constraints.close();

    // allocate memory for [A , B ; B^{T} , 0 ]
    BlockDynamicSparsityPattern dsp(2,2);

    // allocate size for A
    dsp.block(0,0).reinit(n_vector_field,
                          n_vector_field);

    //Allocate size for B^{T}
    dsp.block(1,0).reinit (n_potential,
                           n_vector_field);

    // allocate size for B
    dsp.block(0,1).reinit (n_vector_field,
                           n_potential);

    // allocate size for 0
    dsp.block(1,1).reinit (n_potential,
                           n_potential);

    dsp.collect_sizes();

    // creat the actual sparsity pattern:
    DoFTools::make_sparsity_pattern (dof_handler,
                                     dsp,
                                     constraints);

    constraints.condense(dsp);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit (sparsity_pattern);

    // allocate memory for solutions
    solution.reinit (2); // [vector field, potential]
    solution.block(0).reinit (n_vector_field); // VECTOR FIELD
    solution.block(1).reinit (n_potential); // POTENTIAL

    solution.collect_sizes ();

    // memeory for RHS
    system_rhs.reinit (2);
    system_rhs.block(0).reinit (n_vector_field); // DIRICHLET BC
    system_rhs.block(1).reinit (n_potential); // RIGHT HAND SIDE
    system_rhs.collect_sizes ();

}

template <int dim>
void MixedPoissonProblem<dim>::enforce_Neumann_boundary_values()
{
    const FEValuesExtractors::Vector VectorField(0); // Vector as in Vector field
    ComponentMask vector_field_mask =	fe.component_mask(VectorField);

    DoFTools::make_zero_boundary_constraints(dof_handler,
            1,
            constraints,
            vector_field_mask);
    constraints.close();
}

template<int dim>
void
MixedPoissonProblem<dim>::
assemble_local_matrix(const typename DoFHandler<dim>::active_cell_iterator & cell,
                      Assembly::AssemblyScratch<dim>						& scratch,
                      Assembly::CopyData<dim>										& data,
                      const double 															& lambda)
{
    const unsigned int	dofs_per_cell 	=	scratch.fe_values.dofs_per_cell;
    const unsigned int 	n_q_points 			= scratch.fe_values.n_quadrature_points;

    cell->get_dof_indices(data.local_dof_indices);


    // Get the actual values for vector field and potential from FEValues
    // Use Extractors instead of having to deal with shapefunctions directly
    const FEValuesExtractors::Vector VectorField(0); // Vector as in Vector field
    const FEValuesExtractors::Scalar Potential(dim);

    // get fe_values for this cell
    scratch.fe_values.reinit(cell);
    data.local_matrix = 0;

    // BODY INTEGRALS
    for(unsigned int q=0; q<n_q_points; q++)
    {
        //  Test function basis
        for(unsigned int i=0; i<dofs_per_cell; i++)
        {
            // get the tensor values of the i-th VectorField basis functions
            // at point q
            const Tensor<1, dim>	psi_i_field = scratch.fe_values[VectorField].value(i,q);

            // get the value of the div. of the i-th VectorField basis functions
            // at the point q
            const double	div_psi_i_field =	scratch.fe_values[VectorField].divergence(i,q);

            // get the value of the i-th potential basis functions at the point q
            const double psi_i_potential = scratch.fe_values[Potential].value(i,q);

            // trial functions basis
            for(unsigned int j=0; j<dofs_per_cell; j++)
            {
                // get the tensor values of the j-th VectorField basis functions
                // at point q
                const Tensor<1, dim>	psi_j_field = scratch.fe_values[VectorField].value(j,q);

                // get the value of the div. of the j-th VectorField basis functions
                // at the point q
                const double	div_psi_j_field =	scratch.fe_values[VectorField].divergence(j,q);

                // get the value of the i-th potential basis functions at the point q
                const double psi_j_potential = scratch.fe_values[Potential].value(j,q);

                // build whole matrix at once, not blocks individually.
                // \int (P * D - div P * phi - v * D ) dx
                data.local_matrix(i,j) += ((1.0/lambda) * psi_i_field
                                           * psi_j_field
                                           - div_psi_i_field * psi_j_potential
                                           - psi_i_potential * div_psi_j_field
                                          ) * scratch.fe_values.JxW(q);

            } // for j
        } // for i
    } // for q
} // assemble_local_matrix


template<int dim>
void
MixedPoissonProblem<dim>::
assemble_local_rhs(const typename DoFHandler<dim>::active_cell_iterator & cell,
                   Assembly::AssemblyScratch<dim>						& scratch,
                   Assembly::CopyData<dim>								& data)
{
    const unsigned int	dofs_per_cell 	=	scratch.fe_values.dofs_per_cell;
    const unsigned int 	n_q_points 			= scratch.fe_values.n_quadrature_points;
    const unsigned int	n_face_q_points =	scratch.fe_face_values.n_quadrature_points;

    cell->get_dof_indices(data.local_dof_indices);

    const RightHandSide<dim>							right_hand_side;
    const DirichletBoundaryValues<dim>  	dirichlet_boundary_values;

    // Get the actual values for vector field and potential from FEValues
    // Use Extractors instead of having to deal with shapefunctions directly
    const FEValuesExtractors::Vector VectorField(0); // Vector as in Vector field
    const FEValuesExtractors::Scalar Potential(dim);

    // get fe_values for this cell
    scratch.fe_values.reinit(cell);
    data.local_rhs =	0;

    // get RHS values on this cell
    right_hand_side.value_list(scratch.fe_values.get_quadrature_points(),
                               scratch.rhs_values);

    // BODY INTEGRALS
    for(unsigned int q=0; q<n_q_points; q++)
    {
        //  Test function basis
        for(unsigned int i=0; i<dofs_per_cell; i++)
        {

            // get the value of the i-th potential basis functions at the point q
            const double psi_i_potential = scratch.fe_values[Potential].value(i,q);

            // get the local RHS values
            data.local_rhs(i) += -psi_i_potential *
                                 scratch.rhs_values[q] *
                                 scratch.fe_values.JxW(q);
        } // for i
    } // for q


    // BOUNDARY INTEGRALS
    for(unsigned int face_no=0;
            face_no < GeometryInfo<dim>::faces_per_cell;
            face_no++)
    {
        typename DoFHandler<dim>::face_iterator 	face = cell->face(face_no);

        // make sure only on Dirichlet segmant
        if((face->at_boundary()) && (face->boundary_id() == 0) )
        {
            // get the values of the shape functions at this boundary face
            scratch.fe_face_values.reinit(cell,face_no);

            dirichlet_boundary_values.value_list(
                scratch.fe_face_values.get_quadrature_points(),
                scratch.bc_values);

            for(unsigned int q=0; q<n_face_q_points; q++)
            {
                for(unsigned int i=0; i<dofs_per_cell; i++)
                {
                    // get the i-th VectorFields basis function value at point q
                    data.local_rhs(i)	+=	-(scratch.fe_face_values[VectorField].value(i,q) *
                                              scratch.fe_face_values.normal_vector(q) *
                                              scratch.bc_values[q] *
                                              scratch.fe_face_values.JxW(q));
                } // for i
            } // for q
        } // end if
    } // end for face_no
} // assemble_local_rhs

template<int dim>
void
MixedPoissonProblem<dim>::
copy_local_to_global_matrix(const Assembly::CopyData<dim> & data)
{
    constraints.distribute_local_to_global(data.local_matrix,
                                           data.local_dof_indices,
                                           system_matrix);

} // copy_local_to_global_matrix

template<int dim>
void
MixedPoissonProblem<dim>::
copy_local_to_global_rhs(const Assembly::CopyData<dim> & data)
{
    constraints.distribute_local_to_global(data.local_rhs,
                                           data.local_dof_indices,
                                           system_rhs);


} // copy_local_to_global_matrix


template <int dim>
void MixedPoissonProblem<dim>::serial_assemble_matrix()
{
    Assembly::AssemblyScratch<dim> 	scratch(fe,
                                            QGauss<dim>(degree+2), // body quad
                                            QGauss<dim-1>(degree+2)); // BC quad

    Assembly::CopyData<dim>				data(fe);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    const double lambda = 1.0;

    auto assemble_bound_local_matrix =
        std_cxx11::bind(&MixedPoissonProblem::assemble_local_matrix,
                        this,
                        std_cxx11::_1,
                        std_cxx11::_2,
                        std_cxx11::_3,
                        lambda);

    // LOOP OVER ALL THE ACTIVE CELLS AND BUILD THE BLOCK MATRIX AND RHS
    for(; cell != endc ; cell++)
    {
        assemble_bound_local_matrix(cell, scratch, data);
        copy_local_to_global_matrix(data);
    } // cell

} // assemble_matrix



template <int dim>
void MixedPoissonProblem<dim>::serial_assemble_rhs()
{
    Assembly::AssemblyScratch<dim> 	scratch(fe,
                                            QGauss<dim>(degree+2), // body quad
                                            QGauss<dim-1>(degree+2)); // BC quad

    Assembly::CopyData<dim> 				data(fe);


    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    // LOOP OVER ALL THE ACTIVE CELLS AND BUILD THE BLOCK MATRIX AND RHS
    for(; cell != endc ; cell++)
    {
        assemble_local_rhs(cell, scratch, data);
        copy_local_to_global_rhs(data);

    } // end for cell

} // end assemble_system()



template<int dim>
void
MixedPoissonProblem<dim>::
parallel_assemble_matrix()
{

    const double lambda = 1.0;

    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    std_cxx11::bind(&MixedPoissonProblem<dim>::assemble_local_matrix,
                                    this,
                                    std_cxx11::_1,
                                    std_cxx11::_2,
                                    std_cxx11::_3,
                                    lambda),
                    std_cxx11::bind(&MixedPoissonProblem<dim>::copy_local_to_global_matrix,
                                    this,
                                    std_cxx11::_1),
                    Assembly::AssemblyScratch<dim>(fe,
                            QGauss<dim>(degree+2),
                            QGauss<dim-1>(degree+2)),
                    Assembly::CopyData<dim>(fe)
                   );


} // parallel_assemble_matrix


template<int dim>
void
MixedPoissonProblem<dim>::
parallel_assemble_rhs()
{
    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    *this,
                    &MixedPoissonProblem::assemble_local_rhs,
                    &MixedPoissonProblem::copy_local_to_global_rhs,
                    Assembly::AssemblyScratch<dim>(fe,
                            QGauss<dim>(degree+2),
                            QGauss<dim-1>(degree+2)),
                    Assembly::CopyData<dim>(fe)
                   );
}

template <int dim>
void MixedPoissonProblem<dim>::solve()
{
    direct_solver.initialize(system_matrix);
    direct_solver.vmult(solution,system_rhs);
    constraints.distribute(solution);

} // solve()

template <int dim>
void MixedPoissonProblem<dim>::output_results() const
{
    std::vector<std::string> solution_names;
    switch(dim)
    {
    case 2:
        solution_names.push_back("E_x");
        solution_names.push_back("E_y");
        solution_names.push_back("Pot");
        break;

    case 3:
        solution_names.push_back("E_x");
        solution_names.push_back("E_y");
        solution_names.push_back("E_z");
        solution_names.push_back("Pot");
        break;

    default:
        Assert(false, ExcNotImplemented() );
    }

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, solution_names);

    data_out.build_patches(degree+1);
    std::ofstream output("solution.vtk");
    data_out.write_vtk(output);
}

template <int dim>
void MixedPoissonProblem<dim>::output_system() const
{
    std::ofstream output_system("A.mtx");
    system_matrix.print_formatted(output_system);
    output_system.close();
    output_system.open("b.vec");
    system_rhs.print(output_system);
    output_system.close();
}

template <int dim>
void MixedPoissonProblem<dim>::compute_errors(double & pot_l2_error,
        double & vec_field_l2_error)
{
    const ComponentSelectFunction<dim> potential_mask(dim, dim+1);
    const ComponentSelectFunction<dim>
    vectorField_mask(std::make_pair(0,dim), dim+1);

    TrueSolution<dim>	exact_solution;
    Vector<double> cellwise_errors(triangulation.n_active_cells() );

    QTrapez<1>				q_trapez;
    QIterated<dim> 		quadrature(q_trapez, degree+2);


    VectorTools::integrate_difference(dof_handler, solution, exact_solution,
                                      cellwise_errors, quadrature,
                                      VectorTools::L2_norm,
                                      &potential_mask);

    pot_l2_error = cellwise_errors.l2_norm();

    VectorTools::integrate_difference(dof_handler, solution, exact_solution,
                                      cellwise_errors, quadrature,
                                      VectorTools::L2_norm,
                                      &vectorField_mask);

    vec_field_l2_error = cellwise_errors.l2_norm();

    std::cout << "\nErrors: ||e_pot||_L2 = " << pot_l2_error
              << ",		||e_vec_field||_L2 = " << vec_field_l2_error
              << std::endl;

}

template <int dim>
void MixedPoissonProblem<dim>::run()
{
    /*		make_grid_and_dofs();
    		make_Neumann_boundaries();
    		enforce_Neumann_boundary_values();
    		serial_assemble_matrix();
    		serial_assemble_rhs();
    //		output_system();
    		solve();
    		output_results();
    */
}

template <int dim>
void MixedPoissonProblem<dim>::run_with_errors(double & h,
        double & pot_errors,
        double & vec_field_errors)
{
    TimerOutput timer(std::cout, TimerOutput::summary,
                      TimerOutput::wall_times);

    timer.enter_subsection("Boiler Plate");
    Grid<dim> grid;
    grid.make_local_refined_grid(triangulation,
                                 global_refinement_number,
                                 local_refinement_number);
    grid.make_Dirichlet_boundaries(triangulation);
    grid.print_grid(triangulation);
    make_dofs();
    //make_Neumann_boundaries();
//		enforce_Neumann_boundary_values();
    timer.leave_subsection("Boiler Plate");


    timer.enter_subsection("Assemble Matrix");
//		serial_assemble_matrix();
    parallel_assemble_matrix();
    timer.leave_subsection("Assemble Matrix");

    timer.enter_subsection("Assemble RHS");
//		serial_assemble_rhs();
    parallel_assemble_rhs();
    timer.leave_subsection("Assemble RHS");

//		output_system();

    timer.enter_subsection("Solve");
    solve();
    timer.leave_subsection("Solve");

    timer.enter_subsection("Printing");
    output_results();
    timer.leave_subsection("Printing");

    compute_errors(pot_errors, vec_field_errors);
    h = GridTools::maximal_cell_diameter(triangulation);

}

} // NAMESPACE MIXEDFEM



///////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION
///////////////////////////////////////////////////////////////////////////////


