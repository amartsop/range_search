#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>

#include "geom.h"
#include "geometry2D_mesh.h"
#include "gq_triangle_2D.hpp"
#include "gq_line_1D.hpp"

class Model 
{
    public:
        Model() {};

        // Set integration mesh
        void set_integration_mesh(const Geometry2DMesh::Triangular2DMesh& mesh);

        // Integrate function
        void integrate(void);

        // Integrand function
        double f(const arma::dvec& x);

    private:

        // Gauss quadrature triangle handle
        GQTriangle2D<Model> m_gq_triangle;

        // Gauss quadrature line handle
        GQLine1D<Model> m_gq_line;

        // Integration mesh
        Geometry2DMesh::Triangular2DMesh m_integr_mesh;
};