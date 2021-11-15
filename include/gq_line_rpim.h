#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>

#include "./geom.h"
#include "./geometry2D_mesh.h"
#include "./mesh2D.h"

class GQLineRPIM
{

public:
    // Constructor
    GQLineRPIM();

    // Generate quadrature points 
    geom::PointCloud<double> generate_quadrature_points(const Mesh2D& mesh,
        int evals=4);

    // Quadrature points per element
    struct CellProperties
    {
        // Element idx
        size_t cell_idx;

        // Element length
        double length;

        // Quadrature points indices
        std::vector<size_t> quadr_pt_idx;
    };

    // Get cell properties
    std::vector<CellProperties> get_integration_cells_properties(void);

    // Get quadrature weights
    std::vector<double> get_quadrature_weights(void) { return m_weights; }

private:
    // Weights for given order
    std::vector<double> m_weights;

    // Weights container
    std::vector<std::vector<double>> m_weights_container;

    // Points container
    std::vector<std::vector<double>> m_points_container;

    // Weights and points generation
    void weight_point_generation(void);

    //Quadrature points (x elements - y quadrature points)
    std::vector<std::vector<arma::dvec>> m_quadrature_pts;

    // Elements number
    size_t m_elements_num;

private:
    std::vector<CellProperties> m_cell_properties;
};