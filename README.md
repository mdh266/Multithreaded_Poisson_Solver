# Multithreaded Poisson Equation Solver

## Introduction
This code is designed to numerically solve the <a href="https://en.wikipedia.org/wiki/Poisson's_equation"> Poisson equation</a> using the <a href="https://en.wikipedia.org/wiki/Mixed_finite_element_method"> mixed finite element method</a>.  The code runs in parallel using multithreading through the Intel Thread Building Blocks.

**Note** 
This project improves upon <a href="https://www.dealii.org/8.4.1/doxygen/deal.II/step_20.html"> step-20 </a>in the deal.ii tutorial by:

- Adding Neumann boundary conditions.
- Allow for multithreading to reduce runtimes.
 
## Requirements
The requirements for this software is deal.ii library version 8.4.0 or highe and CMake version 2.8 or higher.

## Installation
First obtain and install a copy of the dealii deal.ii library version 8.4.0 or higher. 

## Compiling
To generate a makefile for this code using CMake type into the terminal:

*cmake . -DDEAL_II_DIR=/path_to_deal.ii*

To compile the code in release mode use:

*make release*

## Running
To run the executable use:

*./main*
