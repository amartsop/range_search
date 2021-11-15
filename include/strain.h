#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>
#include "shape_function.h"

class Strain
{

public:
    Strain() {};
    
    // Set shape function
    void set_shape_function(const ShapeFunction::Measures& sf);

    // Update function
    void update(const arma::dvec& es);

    // Get strain vector
    arma::dvec get_strain_vector(void) { return m_strain_vector; }

    // Get Ds matrix
    arma::dmat get_ds_matrix(void) { return m_ds_matrix; }

private:
    // Shape function s handle
    ShapeFunction::Measures m_sf_s;

    // Number of support domain points
    size_t m_ns;

private:
    // Strain vector
    arma::dvec m_strain_vector;

    // Ds matrix 
    arma::dmat m_ds_matrix;
};



