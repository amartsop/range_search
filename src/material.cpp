#include "../include/material.h"

arma::dmat Material::get_elasticity_matrix(void)
{
    // Benchmark material
    // Poisson ratio
    double gamma = 0.3;

    // Young's modulus (N / m^2)
    double epsilon = 3.0e7; 

    // Initialize matrix
    arma::dmat d_matrix = arma::zeros<arma::dmat>(3, 3);
    d_matrix(0, 0) = 1.0; d_matrix(0, 1) = gamma;
    d_matrix(1, 0) = gamma; d_matrix(1, 1) = 1.0;
    d_matrix(2, 2) = (1.0 - gamma) / 2.0;
    d_matrix *= epsilon / (1.0 - pow(gamma, 2.0));

    return d_matrix;
}