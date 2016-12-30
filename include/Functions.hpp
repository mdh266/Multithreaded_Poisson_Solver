#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
namespace MixedFEM
{
using namespace dealii;

template <int dim>
class RightHandSide : public Function<dim>
{
public:
    RightHandSide() : Function<dim>(1)
    {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0 ) const;
};

template <int dim>
class DirichletBoundaryValues : public Function<dim>
{
public:
    DirichletBoundaryValues() : Function<dim>(1)
    {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0 ) const;
};


// Inverse Lambda^{2} Tensor
template<int dim>
class DebeyeInverse : public TensorFunction<2,dim>
{
public:
    DebeyeInverse() : TensorFunction<2,dim>()
    {}

    virtual void value_list(const std::vector<Point<dim> > &points,
                            std::vector<Tensor<2,dim> > &values) const;
};

///////////////////////////////////////////////////////////////////////////////
// TEST CASE
///////////////////////////////////////////////////////////////////////////////
template<int dim>
class TrueSolution : public Function<dim>
{
public:
    TrueSolution() : Function<dim>(dim+1)
    {}

    virtual void vector_value(const Point<dim> & p,
                              Vector<double> &valuess) const;
};

}
#endif
