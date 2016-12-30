#include "../include/Grid.hpp"


/////////////////////////////////////////////////////////////////////////////////
// GRID AND BOUNDARIES
/////////////////////////////////////////////////////////////////////////////////

using namespace dealii;

template<int dim>
Grid<dim>::
Grid()
{
}

template <int dim>
void
Grid<dim>::
make_grid(Triangulation<dim> & triangulation,
          const unsigned int & global_refine)

{
    // make the triangulation and refine globally n_refine times
    GridGenerator::hyper_cube(triangulation,0,1);
    triangulation.refine_global(global_refine);


}

template <int dim>
void
Grid<dim>::
make_test_grid(Triangulation<dim>           & triangulation,
               const int                   & n_global_refine)
{
    // make the triangulation and refine globally n_refine times
    GridGenerator::hyper_cube(triangulation,0,1);
    triangulation.refine_global(n_global_refine);

    // set the boundaries to be Dirichlet
    typename Triangulation<dim>::cell_iterator
    cell = triangulation.begin(),
    endc = triangulation.end();

    // loop over all the cells
    for(; cell != endc; cell++)
    {
        // loop over all the faces of the cell and find which are on the boundary
        for(unsigned int face_no=0;
                face_no < GeometryInfo<dim>::faces_per_cell;
                face_no++)
        {
            if(cell->face(face_no)->at_boundary() )
            {
                if((cell->face(face_no)->center()[1] == 0) ||
                        (cell->face(face_no)->center()[1] == 1.0) )
                {
                    // set it to be Neumann boundary condition by setting the boundary
                    // indicator to be 1.  NOTE: Default is 0, which implies Dirichlet
                    cell->face(face_no)->set_boundary_id(Neumann);
                }
                else
                {
                    //NOTE: Default is 0, which implies Interface
                    cell->face(face_no)->set_boundary_id(Dirichlet);
                }
            } // end if on boundary
        } // for face_no
    } // for cell
} // make test_grid


template <int dim>
void
Grid<dim>::
make_local_refined_grid(Triangulation<dim> & triangulation,
                        const unsigned int & global_refine,
                        const unsigned int & local_refine)
{
    // make the triangulation and refine globally n_refine times
    GridGenerator::hyper_cube(triangulation,0,1);
    triangulation.refine_global(global_refine);

    for(unsigned int i =0; i <local_refine; i++)
    {
        typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active(),
        endc = triangulation.end();

        // loop over cells and locally mark appropriate cells for refinement
        for(; cell != endc; cell++)
        {
//          if(cell->at_boundary())
            if((cell->center()[1]) > 0.9 )
            {
                if((cell->center()[0] > 0.9)  || (cell->center()[0] < 0.1) )
                    cell->set_refine_flag();
            }
        }
        // refine the marked cells
        triangulation.execute_coarsening_and_refinement();
    }


}

template<int dim>
void
Grid<dim>::
make_Dirichlet_boundaries(Triangulation<dim> & triangulation)
{
    // NOTE::   Sets the boundary to be default Dirichlet
    typename Triangulation<dim>::cell_iterator
    cell = triangulation.begin(),
    endc = triangulation.end();

    // loop over all the cells
    for(; cell != endc; cell++)
    {
        // loop over all the faces of the cell and find which are on the boundary
        for(unsigned int face_no=0;
                face_no < GeometryInfo<dim>::faces_per_cell;
                face_no++)
        {
            if(cell->face(face_no)->at_boundary() )
            {
                //NOTE: Default is 0, which implies Dirichlet
                cell->face(face_no)->set_boundary_id(Dirichlet);
            } // end if on boundary
        } // for face_no
    } // for cell

}


template<int dim>
void
Grid<dim>::
make_Neumann_boundaries(Triangulation<dim> & triangulation)
{
    // NOTE:: BECAREFUL TO MATCH WITH CUBE END POINTS
    //        BECAREFUL WITH DIMENSION

    typename Triangulation<dim>::cell_iterator
    cell = triangulation.begin(),
    endc = triangulation.end();
    // loop over all the cells
    for(; cell != endc; cell++)
    {
        // loop over all the faces of the cell and find which are on the bondary
        for(unsigned int face_no=0;
                face_no < GeometryInfo<dim>::faces_per_cell;
                face_no++)
        {
            if(cell->face(face_no)->at_boundary() )
            {
                // set the portions of the boundary y = 0 and y = +1
                // to be Neumann boundary conditions
                if( (cell->face(face_no)->center()[0] == 0) ||
                        (cell->face(face_no)->center()[0] == +1) )
                {
                    // set it to be Neumann boundary condition by setting the boundary
                    // indicator to be 1.  NOTE: Default is 0, which implies Dirichlet
                    cell->face(face_no)->set_boundary_id(Neumann);
                } // end if Neumann segment
            } // end if on boundary
        } // for face_no
    } // for cell

}

template<int dim>
void
Grid<dim>::
print_grid(Triangulation<dim> & triangulation)
{
    std::ofstream out("grid.eps");
    GridOut grid_out;
    grid_out.write_eps(triangulation,out);
    out.close();
}



