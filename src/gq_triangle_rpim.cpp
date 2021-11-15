#include "../include/gq_triangle_rpim.h"

// Constructor
GQTriangleRPIM::GQTriangleRPIM()
{
    // Initialize gaussian points and weights
    weight_point_generation();
}


// Generate quadrature points 
geom::PointCloud<double> GQTriangleRPIM::generate_quadrature_points(const
    Mesh2D& mesh, int order)
{
    // Get weights and points for the order provided
    m_wps = m_wp.at(order);

    // Get quadrature weights
    for (size_t i = 0; i < m_wps.weight.size(); i++)
    {
        m_weights.push_back(m_wps.weight.at(i).wi);
    }

    // Get elements number
    m_elements_num = mesh.volume_elements.size();

    // Initialize quadrature points container
    geom::PointCloud<double> quadr_pts;

    // Initialize quadrature point counter
    size_t quadr_counter = 0;

    // Calculate integral over surface area
    for (size_t i = 0; i < m_elements_num ; i++)
    {
        // Get the vertex coordinates of triangle i
        int index1 = mesh.volume_elements.at(i).node1_index;
        int index2 = mesh.volume_elements.at(i).node2_index;
        int index3 = mesh.volume_elements.at(i).node3_index;

        geom::Point<double> vertex1_i = mesh.node_coords.pts.at(index1);
        geom::Point<double> vertex2_i = mesh.node_coords.pts.at(index2);
        geom::Point<double> vertex3_i = mesh.node_coords.pts.at(index3);

        double x1_i = vertex1_i.x; double y1_i = vertex1_i.y; 
        double x2_i = vertex2_i.x; double y2_i = vertex2_i.y; 
        double x3_i = vertex3_i.x; double y3_i = vertex3_i.y; 

        // Define x3_i_vec
        arma::dvec x3_i_vec = {x3_i, y3_i};

        // Initialize quadrature point indices
        CellProperties cell_props;

        // Set element id
        cell_props.cell_idx = i;

        // Calculate the area of triangle i
        arma::dmat ai_mat = {{1.0, x1_i, y1_i}, {1.0, x2_i, y2_i}, {1.0, x3_i, y3_i}};
        double ai = (1.0 / 2.0) * std::abs(arma::det(ai_mat));
        cell_props.area = ai;
        
        // Jacobian matrix
        arma::dmat ju_i = {{(x1_i - x3_i), (x2_i - x3_i)},
            {(y1_i - y3_i), (y2_i - y3_i)}};

        // Initialize internal vector
        std::vector<arma::dvec> xk_vec;

        for (size_t k = 0; k < m_wps.weight.size(); k++)
        {
            // Get integration points 
            arma::dvec uk = {m_wps.point.at(k).u, m_wps.point.at(k).v};

            // Map integration points
            arma::dvec xk_i_vec = x3_i_vec + ju_i * uk;

            // Push xk_i_vec to quadrature points
            xk_vec.push_back(xk_i_vec);

            // Generate quadrature points point cloud
            geom::Point<double> point_i{xk_i_vec(0), xk_i_vec(1), 0.0};
            quadr_pts.pts.push_back(point_i);

            // Push back quadrature index
            cell_props.quadr_pt_idx.push_back(quadr_counter);

            // Update quadrature point counter
            quadr_counter++;
        }

        // Push quadrature indices to structure
        m_cell_properties.push_back(cell_props);

        // Push xk_vec to quadrature points container
        m_quadrature_pts.push_back(xk_vec);
    }
    
    return quadr_pts;
}

// Get quadrature indices structure
std::vector<GQTriangleRPIM::CellProperties>
    GQTriangleRPIM::get_integration_cells_properties(void)
{
    return m_cell_properties;
}

// Generation of weights and points for gaussian quadrature over 2D triangles
// (Zienkiewicz - The finite element method Volume 1 - Table 9.2)
void  GQTriangleRPIM::weight_point_generation(void)
{
    /************************** Order: Linear **************************/
    // Create weight-point struct and initialize internal vectors
    WeightPoint wp_linear;

    wp_linear.weight.resize(1);
    wp_linear.point.resize(1);

    // Point 1
    wp_linear.weight.at(0).wi = 1.0;
    wp_linear.point.at(0).u = 1.0 / 3.0;
    wp_linear.point.at(0).v = 1.0 / 3.0;
    wp_linear.point.at(0).w = 1.0 / 3.0;

    // Push to vector
    m_wp.push_back(wp_linear);

    /************************** Order: Quadratic **************************/
    // Create weight-point struct and initialize internal vectors
    WeightPoint wp_quadr;

    wp_quadr.weight.resize(3);
    wp_quadr.point.resize(3);

    // Point 1
    wp_quadr.weight.at(0).wi = 1.0 / 3.0;
    wp_quadr.point.at(0).u = 1.0 / 2.0;
    wp_quadr.point.at(0).v = 1.0 / 2.0;
    wp_quadr.point.at(0).w = 0.0;

    // Point 2
    wp_quadr.weight.at(1).wi = 1.0 / 3.0;
    wp_quadr.point.at(1).u = 0.0;
    wp_quadr.point.at(1).v = 1.0 / 2.0;
    wp_quadr.point.at(1).w = 1.0 / 2.0;

    // Point 3
    wp_quadr.weight.at(2).wi = 1.0 / 3.0;
    wp_quadr.point.at(2).u = 1.0 / 2.0;
    wp_quadr.point.at(2).v = 0.0;
    wp_quadr.point.at(2).w = 1.0 / 2.0;

    // Push to vector
    m_wp.push_back(wp_quadr);

    /************************** Order: Cubic **************************/
    // Create weight-point struct and initialize internal vectors
    WeightPoint wp_cub;

    wp_cub.weight.resize(4);
    wp_cub.point.resize(4);

    // Point 1
    wp_cub.weight.at(0).wi = - 27.0 / 48.0;
    wp_cub.point.at(0).u = 1.0 / 3.0;
    wp_cub.point.at(0).v = 1.0 / 3.0;
    wp_cub.point.at(0).w = 1.0 / 3.0;

    // Point 2
    wp_cub.weight.at(1).wi = 25.0 / 48.0;
    wp_cub.point.at(1).u = 0.6;
    wp_cub.point.at(1).v = 0.2;
    wp_cub.point.at(1).w = 0.2;

    // Point 3
    wp_cub.weight.at(2).wi = 25.0 / 48.0;
    wp_cub.point.at(2).u = 0.2;
    wp_cub.point.at(2).v = 0.6;
    wp_cub.point.at(2).w = 0.2;

    // Point 4
    wp_cub.weight.at(3).wi = 25.0 / 48.0;
    wp_cub.point.at(3).u = 0.2;
    wp_cub.point.at(3).v = 0.2;
    wp_cub.point.at(3).w = 0.6;

    // Push to vector
    m_wp.push_back(wp_cub);

    /************************** Order: Quintic **************************/
    // Create weight-point struct and initialize internal vectors
    WeightPoint wp_quint;

    // Define constants
    double a1 = 0.0597158717; 
    double b1 = 0.4701420641; 
    double a2 = 0.797426985; 
    double b2 = 0.1012865073;

    double w1 = 0.1323941527;
    double w2 = 0.1259391805;

    wp_quint.weight.resize(7);
    wp_quint.point.resize(7);

    // Point 1
    wp_quint.weight.at(0).wi = 0.225;
    wp_quint.point.at(0).u = 1.0 / 3.0;
    wp_quint.point.at(0).v = 1.0 / 3.0;
    wp_quint.point.at(0).w = 1.0 / 3.0;

    // Point 2
    wp_quint.weight.at(1).wi = w1;
    wp_quint.point.at(1).u = a1;
    wp_quint.point.at(1).v = b1;
    wp_quint.point.at(1).w = b1;

    // Point 3
    wp_quint.weight.at(2).wi = w1;
    wp_quint.point.at(2).u = b1;
    wp_quint.point.at(2).v = a1;
    wp_quint.point.at(2).w = b1;

    // Point 4
    wp_quint.weight.at(3).wi = w1;
    wp_quint.point.at(3).u = b1;
    wp_quint.point.at(3).v = b1;
    wp_quint.point.at(3).w = a1;

    // Point 5
    wp_quint.weight.at(4).wi = w2;
    wp_quint.point.at(4).u = a2;
    wp_quint.point.at(4).v = b2;
    wp_quint.point.at(4).w = b2;

    // Point 6
    wp_quint.weight.at(5).wi = w2;
    wp_quint.point.at(5).u = b2;
    wp_quint.point.at(5).v = a2;
    wp_quint.point.at(5).w = b2;

    // Point 7
    wp_quint.weight.at(6).wi = w2;
    wp_quint.point.at(6).u = b2;
    wp_quint.point.at(6).v = b2;
    wp_quint.point.at(6).w = a2;

    // Push to vector
    m_wp.push_back(wp_quint);
}