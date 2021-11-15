#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>

#include "./geom.h"
#include "./mesh2D.h"
#include "./geometry2D_mesh.h"

class GQTriangleRPIM
{

public:
    // Integration order
    inline static const int LINEAR = 0;
    inline static const int QUADRATIC = 1;
    inline static const int CUBIC = 2;
    inline static const int QUINTIC = 3;

public:
    // Constructor
    GQTriangleRPIM();

    // Generate quadrature points 
    geom::PointCloud<double> generate_quadrature_points(const
        Mesh2D& mesh, int order = LINEAR);

    // Quadrature points per element
    struct CellProperties
    {
        // Element idx
        size_t cell_idx;

        // Element Area
        double area;

        // Quadrature points indices
        std::vector<size_t> quadr_pt_idx;
    };

    // Get cell properties
    std::vector<CellProperties> get_integration_cells_properties(void);

    // Get quadrature weights
    std::vector<double> get_quadrature_weights(void) { return m_weights;}

private:

    // Struct of weight (Gaussian quadrature)
    struct Weight{ double wi; };

    // Struct of normalised point
    struct Point { double u, v, w; };

    // Weight-point combination
    struct WeightPoint {
        std::vector<Weight> weight;
        std::vector<Point> point; };

    // Define weight-points handle
    std::vector<WeightPoint> m_wp;

    // Weights and points generation
    void weight_point_generation(void);

    //Quadrature points (x elements - y quadrature points)
    std::vector<std::vector<arma::dvec>> m_quadrature_pts;

    // Weighs and points
    WeightPoint m_wps;

    // Elements number
    size_t m_elements_num;

    // Quadrature weights
    std::vector<double> m_weights;

private:
    std::vector<CellProperties> m_cell_properties;
};