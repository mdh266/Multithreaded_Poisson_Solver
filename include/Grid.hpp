#ifndef _GRID_H__
#define _GRID_H__

#include<deal.II/grid/tria.h>
#include<deal.II/grid/tria_accessor.h>
#include<deal.II/grid/tria_iterator.h>
#include<deal.II/grid/grid_generator.h>
#include<deal.II/grid/grid_refinement.h>
#include<deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <iostream>
#include<fstream>


using namespace dealii;


/// \brief This object will be used to build meshes over triangulations provided.
template<int dim>
class Grid
{
public:
    /** \brief Empty constructor.*/
    Grid();

    /** \brief Makes a grid.*/
    void make_grid(Triangulation<dim> & triangulation,
                   const unsigned int & global_refine);

    /**  \brief Makes the grid for the manufactured solution.*/
    void make_test_grid(Triangulation<dim>  & triangulation,
                        const int           & n_global_refine);

    /** \brief Makes the grid for the manufactured solution with locally refined mesh.*/
    void make_local_refined_grid(Triangulation<dim> & triangulation,
                                 const unsigned int & global_refine,
                                 const unsigned int & local_refine);

    /** \brief Loops the cells and marks faces as Dirichlet.
    Faces that are on the Dirchilet boundary are defined here. */
    void make_Dirichlet_boundaries(Triangulation<dim> & triangulation);

    /** \brief Loops the cells and marks faces as Neumann.
    Faces that are on the Neumann boundary are defined here. */
    void make_Neumann_boundaries(Triangulation<dim> & triangulation);

    /** \brief Prints the grid to a "grid.eps". */
    void print_grid(Triangulation<dim> & triangulation);

private:

    /** \brief Enum is used to define the types of boundaries.*/
    enum
    {
        Dirichlet,
        Neumann
    };
};


#endif
