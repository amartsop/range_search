#include "../include/point_loads.h"

void PointLoads::initialize(const Mesh2D& mesh)
{
    // Select nodal coordinates
//    m_selected_nodes_id = node_selection(mesh.node_coords);
    m_selected_nodes_id.push_back(28);
}

// Update point load
arma::dvec PointLoads::get_point_loads(const arma::dvec& q_bar)
{
    // Initialize point loads 
    arma::dvec qp_vector = arma::zeros(q_bar.n_rows, 1);

    // Loop through selected nodes
    for (size_t i = 0; i < m_selected_nodes_id.size(); i++)
    {
        // Get selected node id
        size_t node_i = m_selected_nodes_id.at(i);

        // Get indices of selected node
        size_t index_1 = 2.0 * node_i;
        size_t index_2 = 2.0 * node_i + 1;

        // Calculate point load vector (can change depending on the node)
        arma::dvec point_load = point_load_function(q_bar);

        // Assign point load to the specific nodal coordinates
        qp_vector.at(index_1) = point_load(0);
        qp_vector.at(index_2) = point_load(1);
   }
    return qp_vector;
}

// Point load function
arma::dvec PointLoads::point_load_function(const arma::dvec& q_bar)
{
    // Initialize point load
    arma::dvec point_load = {0.0, 0.0};
    // arma::dvec point_load = {0.0, -1.0e3};
    // arma::dvec point_load = {0.0, -8.0e4};
    // arma::dvec point_load = {0.0, -5.0e5};

    return point_load;
}

// Node selection function
std::vector<int> PointLoads::node_selection(const geom::PointCloud<double>&
    nodes_coords, double tol)
{
    // Initialize target nodes id container
    std::vector<int> node_ids;

    // Get number of nodes
    size_t nodes_num = nodes_coords.pts.size();

    // Convert nodes to vector structure
    auto nodes_vec = geom::conv_pc_to_pc_vec(nodes_coords);
    
    // Find nodes ranges
    auto range_x = gp_utils::find_vector_range(nodes_vec.x);
    auto range_y = gp_utils::find_vector_range(nodes_vec.y);

    // Get target point
    double x_target = range_x.max;
    double y_target = range_y.min + range_y.val / 2.0;
    arma::dvec target_node = {x_target, y_target};

    for (size_t i = 0; i < nodes_num; i++)
    {
        // Create vector of node i
        double node_i_x = nodes_coords.pts.at(i).x;
        double node_i_y = nodes_coords.pts.at(i).y;
        arma::dvec node_i = {node_i_x, node_i_y};

        // Nodal difference
        double nd = arma::norm(target_node - node_i);

        if (nd <= tol)
        {
            node_ids.push_back(i);
        }
    }

    return node_ids;
}