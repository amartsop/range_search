#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>

class ShapeFunction
{

public:
    ShapeFunction(double ac, double dc, double q, int ms=3);

    // Returns the vector and its jacoban
    struct VecJac {
        // Vector
        arma::dvec vec;

        // Jacobian
        arma::dmat jac;
    };

    struct Measures {
        // Phis vector
        arma::dvec phis_vec;

        // Phis jacobian
        arma::dmat phis_jac;

        // Phis Matrix
        arma::dmat phis_mat;
    };

    // Radial basis function (multi-quadrics)
    VecJac rbf_mq(const arma::dvec& x,
        const std::vector<arma::dvec>& sup_dom);

    // Polynomial function (linear basis; m=3)
    VecJac polynomial2D_basis(const arma::dvec& x);
    
    // Calculate 
    Measures calculate(const arma::dvec& x, const std::vector<arma::dvec>& sup_dom);

private:
    // Shape function constants
    double m_ac, m_dc, m_q;

    // Polynomial basis order
    int m_ms;

private: 
    // Gs matrix
    arma::dmat gs_matrix(const std::vector<arma::dvec>& sup_dom);
};


