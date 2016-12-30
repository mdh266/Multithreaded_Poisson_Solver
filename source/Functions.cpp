#include "../include/Functions.hpp"

namespace MixedFEM
{

template <int dim>
void
DebeyeInverse<dim>::value_list(const std::vector<Point<dim> > &points,
                               std::vector<Tensor<2,dim> > &values) const
{
    Assert( points.size() == values.size(),
            ExcDimensionMismatch(points.size(), values.size() ) );

    // returns vector of tensors evaluted at point point[p]
    for(unsigned int p=0; p < points.size(); p++)
    {
        values[p].clear();
        for(unsigned int d=0; d<dim; d++)
        {
            values[p][d][d] = 1.0; // K^{-1}_{dd}(point[p])
        }
    }
}


template <int dim>
double RightHandSide<dim>::value(const Point<dim> &p,
                                 const unsigned int ) const
{
    double x = p[0];
    double z = x*x - 1;
    double y = p[1];

//		double return_value = 0.0;
//		for(unsigned int i=0; i<dim; i++)
    //		return_value += -6.0 * p[i];

    return -6*z*y*(5*x*x-1);
}

template <int dim>
double DirichletBoundaryValues<dim>::value(const Point<dim> &p,
        const unsigned int ) const
{

    /*		double return_value;
    		for(unsigned int d=0; d<dim; d++)
    			return_value += p[d]*p[d]*p[d];

    		return return_value;
    */
    double x = p[0];
    double z = x*x - 1;
    double y = p[1];

    return z*z*z*y;

}



template <int dim>
void TrueSolution<dim>::vector_value(const Point<dim> &p,
                                     Vector<double> &values) const
{
    Assert(values.size() == dim+1,
           ExcDimensionMismatch(values.size(), dim+1) );

    double x = p[0];
    double y = p[1];
    double z = (x*x - 1);

    values(0) = -6*x*z*z*y;
    values[1] = -z*z*z;
    values[2] =	z*z*z*y;

    /*
    		double pot_value = 0.0;
    		for(unsigned int i = 0; i < dim; i++)
    		{
    			values(i) = -3 *p(i) * p(i);
    			pot_value += p(i) * p(i) * p(i);
    		}
    		values(dim) = pot_value;
    */
}

}