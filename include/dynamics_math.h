#pragma once 
#include <iostream>
#include <armadillo>

namespace dm
{
    struct State{
        // Displacement 
        arma::dvec q;

        // Velocity
        arma::dvec q_dot;

        // Acceleration
        arma::dvec q_ddot;
    };

    // Element id and coordinate
    struct GenForce {
        std::vector<double> Fx; std::vector<double> Fy; std::vector<double> Fz;
        std::vector<double> Mx; std::vector<double> My; std::vector<double> Mz;
    };

    struct Trajectory {
        std::vector<double> x; std::vector<double> y; std::vector<double> z;
        std::vector<double> phi; std::vector<double> theta; std::vector<double> psi;
    };

    // Locator matrix
    arma::dmat locator_matrix(arma::ivec locator_vector, 
        uint32_t n_cols);

    // Skew-symmetric matrix 
    arma::dmat33 s(arma::dvec r);
} 

