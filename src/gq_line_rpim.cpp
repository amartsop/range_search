#include "../include/gq_line_rpim.h"

// Constructor
GQLineRPIM::GQLineRPIM()
{
    // Initialize gaussian points and weights
    weight_point_generation();
}

// Generate quadrature points 
geom::PointCloud<double> GQLineRPIM::generate_quadrature_points(const
    Mesh2D& mesh, int evals)
{
    // Get weights for order provided (evals)
    m_weights = m_weights_container.at(evals - 1);

    // Get elements number
    m_elements_num = mesh.bound_elements.size();

    // Initialize quadrature points container
    geom::PointCloud<double> quadr_pts;

    // Initialize quadrature point counter
    size_t quadr_counter = 0;

    // Calculate integral over line
    for (size_t i = 0; i < m_elements_num; i++)
    {
        // Get the vertex coordinates of line element i
        int index1 = mesh.bound_elements.at(i).node1_index;
        int index2 = mesh.bound_elements.at(i).node2_index;

        geom::Point<double> vertex1_i = mesh.node_coords.pts.at(index1);
        geom::Point<double> vertex2_i = mesh.node_coords.pts.at(index2);

        double x1_i = vertex1_i.x; double y1_i = vertex1_i.y; 
        double x2_i = vertex2_i.x; double y2_i = vertex2_i.y; 

        // Initialize quadrature point indices
        CellProperties cell_props;

        // Set element id
        cell_props.cell_idx = i;

        // Calculate ls
        double ls = (1.0 / 2.0) * std::sqrt(pow(x2_i - x1_i, 2.0) +
            pow(y2_i - y1_i, 2.0));
        cell_props.length = ls;

        // Initialize internal vector
        std::vector<arma::dvec> xk_vec;

        for (size_t k = 0; k < evals; k++)
        {
            double ksi_k = m_points_container.at(evals-1).at(k);

            // Get integration points 
            double xk_i = (x2_i - x1_i) * ksi_k / 2.0 + (x1_i + x2_i) / 2.0;
            double yk_i = (y2_i - y1_i) * ksi_k / 2.0 + (y1_i + y2_i) / 2.0;

            // Map integration points
            arma::dvec xk_i_vec = arma::dvec({xk_i, yk_i});

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
std::vector<GQLineRPIM::CellProperties>
    GQLineRPIM::get_integration_cells_properties(void)
{
    return m_cell_properties;
}

// Generation of weights and points for gaussian quadrature
void  GQLineRPIM::weight_point_generation(void)
{
    // One evaluation
    std::vector<double> weights_1 = {2.0};
    std::vector<double> points_1 = {0.0};
    m_weights_container.push_back(weights_1);
    m_points_container.push_back(points_1);

    // Two evaluations
    std::vector<double> weights_2 = {1.0, 1.0};
    std::vector<double> points_2 = {-0.5773502691896257, 0.5773502691896257};
    m_weights_container.push_back(weights_2);
    m_points_container.push_back(points_2);

    // Three evaluations
    std::vector<double> weights_3 = { 0.8888888888888888, 0.5555555555555556,
        0.5555555555555556};
    std::vector<double> points_3 = {0.0, -0.7745966692414834, 0.7745966692414834};
    m_weights_container.push_back(weights_3);
    m_points_container.push_back(points_3);

    // Four evaluations
    std::vector<double> weights_4 = { 0.6521451548625461, 0.6521451548625461,
        0.3478548451374538, 0.3478548451374538};
    std::vector<double> points_4 = {-0.3399810435848563, 0.3399810435848563,
        -0.8611363115940526, 0.8611363115940526};
    m_weights_container.push_back(weights_4);
    m_weights_container.push_back(points_4);

    // Five evaluations
    std::vector<double> weights_5 = { 0.5688888888888889, 0.4786286704993665,
        0.4786286704993665, 0.2369268850561891, 0.2369268850561891};
    std::vector<double> points_5 = { 0.0, -0.5384693101056831,
        0.5384693101056831, -0.9061798459386640, 0.9061798459386640};
    m_weights_container.push_back(weights_5);
    m_points_container.push_back(points_5);

    // Six evaluations
    std::vector<double> weights_6 = { 0.3607615730481386, 0.3607615730481386,
        0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 
        0.1713244923791704};
    std::vector<double> points_6 = { 0.6612093864662645, -0.6612093864662645,
        -0.2386191860831969, 0.2386191860831969, -0.9324695142031521, 
        0.9324695142031521};
    m_weights_container.push_back(weights_6);
    m_points_container.push_back(points_6);

    // Seven evaluations
    std::vector<double> weights_7 = { 0.4179591836734694, 0.3818300505051189,
        0.3818300505051189, 0.2797053914892766, 0.2797053914892766, 
        0.1294849661688697, 0.1294849661688697};
    std::vector<double> points_7 = { 0.0, 0.4058451513773972,
        -0.4058451513773972, -0.7415311855993945, 0.7415311855993945, 
        -0.9491079123427585, 0.9491079123427585};
    m_weights_container.push_back(weights_7);
    m_points_container.push_back(points_7);

    // Eight evaluations
    std::vector<double> weights_8 = {0.3626837833783620, 0.3626837833783620,
        0.3137066458778873, 0.3137066458778873, 0.2223810344533745,
        0.2223810344533745, 0.1012285362903763, 0.1012285362903763};
    std::vector<double> points_8 = { -0.1834346424956498, 0.1834346424956498,
        -0.5255324099163290, 0.5255324099163290, -0.7966664774136267,
        0.7966664774136267, -0.9602898564975363, 0.9602898564975363};
    m_weights_container.push_back(weights_8);
    m_points_container.push_back(points_8);

    // Nine evaluations
    std::vector<double> weights_9 = { 0.3302393550012598, 0.1806481606948574,
	    0.1806481606948574, 0.0812743883615744, 0.0812743883615744,
        0.3123470770400029, 0.3123470770400029, 0.2606106964029354,
        0.2606106964029354};
    std::vector<double> points_9 = { 0.0000000000000000, -0.8360311073266358,
        0.8360311073266358, -0.9681602395076261, 0.9681602395076261,
        -0.3242534234038089, 0.3242534234038089, -0.6133714327005904,
	    0.6133714327005904};
    m_weights_container.push_back(weights_9);
    m_points_container.push_back(points_9);

    // Ten evaluations
    std::vector<double> weights_10 = {0.2955242247147529, 0.2955242247147529,
        0.2692667193099963, 0.2692667193099963, 0.2190863625159820,
        0.2190863625159820, 0.1494513491505806, 0.1494513491505806,
	    0.0666713443086881, 0.0666713443086881};
    std::vector<double> points_10 = {-0.1488743389816312, 0.1488743389816312,
        -0.4333953941292472, 0.4333953941292472, -0.6794095682990244,
        0.6794095682990244, -0.8650633666889845, 0.8650633666889845,
        -0.9739065285171717, 0.9739065285171717};
    m_weights_container.push_back(weights_10);
    m_points_container.push_back(points_10);
}
