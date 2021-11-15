#include <iostream>
#include <vector>
#include <armadillo>

#include "./geom.h"
#include "./mesh2D.h"
#include "./gp_utils.h"
#include "./dynamics_math.h"

class PointLoads
{

public:
    PointLoads() {}

    // Initilize point loads
    void initialize(const Mesh2D& mesh);

    // Get point load 
    arma::dvec get_point_loads(const arma::dvec& q_bar);

private:

    // Node selection function
    std::vector<int> node_selection(const geom::PointCloud<double>&
        nodes_coords, double tol = 1.0e-6);

    // Selected nodes id
    std::vector<int> m_selected_nodes_id;

    // Point load function
    arma::dvec point_load_function(const arma::dvec& q_bar);
    
};