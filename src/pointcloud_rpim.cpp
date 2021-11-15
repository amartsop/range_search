#include "../include/pointcloud_rpim.h"

void PointcloudRPIM::initialize(const Mesh2D& field_nodes_mesh,
    size_t volume_quadr_interp, size_t surface_quadr_interp)
{
    // Get field nodes mesh
    m_field_nodes_mesh = field_nodes_mesh;

    // Get number of field nodes
    m_field_nodes_num = m_field_nodes_mesh.node_coords.pts.size();

    // Get volume quadrature interpolation number
    m_volume_quadr_interp = volume_quadr_interp;

    // Get surface quadrature interpolation number
    m_surface_quadr_interp = surface_quadr_interp;
    
    // Generate volume quadrature points
    generate_volume_quadrature_points();
    
    // Generate surface quadrature points
    generate_surface_quadrature_points();

    // Get the number of volume quadrature points
    m_volume_gq_pts_num = m_volume_gq_pts.pts.size();

    // Get the number of surface quadrature points
    m_surface_gq_pts_num = m_surface_gq_pts.pts.size();

    // Generate data points structure
    generate_data_points();
}

// Generate data points structure
void PointcloudRPIM::generate_data_points(void)
{
    /****************** Construct global data container *************************/
    // Global points number
    size_t global_pts_num = m_field_nodes_num + m_volume_gq_pts_num +
        m_surface_gq_pts_num;

    // Check quadrature points
    quadrature_points_check();

    // Clear previous values and empty memory
    m_cloud.pts.clear();

    // Reserve values and insert data points
    m_cloud.pts.reserve(global_pts_num);

    // Append field nodes
    m_cloud.pts.insert(m_cloud.pts.end(),
        m_field_nodes_mesh.node_coords.pts.begin(),
        m_field_nodes_mesh.node_coords.pts.end());
        
    // Append volumes's quadrature points
    m_cloud.pts.insert(m_cloud.pts.end(), m_volume_gq_pts.pts.begin(),
        m_volume_gq_pts.pts.end());

    // Append surface's quadrature points
    m_cloud.pts.insert(m_cloud.pts.end(), m_surface_gq_pts.pts.begin(),
        m_surface_gq_pts.pts.end());
}


// Generate triangle quadrature points
void PointcloudRPIM::generate_volume_quadrature_points(void)
{
    // Get quadrature points for volume elements
    m_volume_gq_pts = m_volume_gq.generate_quadrature_points(m_field_nodes_mesh,
        m_volume_quadr_interp);

    // Get volume cell properties
    m_volume_cell_props = m_volume_gq.get_integration_cells_properties();

    // Get volume cell quadrature weights
    m_volume_cell_quadr_weights = m_volume_gq.get_quadrature_weights();
}

// Generate triagle quadrature points
void PointcloudRPIM::generate_surface_quadrature_points(void)
{
    // Get quadrature points for surface elements
    m_surface_gq_pts = m_surface_gq.generate_quadrature_points(m_field_nodes_mesh,
        m_surface_quadr_interp);

    // Get surface cell properties
    m_surface_cells_props = m_surface_gq.get_integration_cells_properties();

    // Get surface cell quadrature weights
    m_surface_cell_quadr_weights = m_surface_gq.get_quadrature_weights();
}

// Quadrature point check
void PointcloudRPIM::quadrature_points_check(void)
{
    // Total number of field nodes
    double field_nodes_pts = (double) m_field_nodes_num;

    // Total number of quadrature points
    double quadrature_pts = (double) m_volume_gq_pts_num +
        (double) m_surface_gq_pts_num;

    // Quadrature to field nodes ratio
    double ratio = quadrature_pts / field_nodes_pts;

    if(ratio >= 3.0)
    {
        std::cout << "Quadrature points number: " << std::to_string(ratio) << ", OK!" << std::endl;
    }
    else
    {
        std::cout << "Quadrature points number: " << std::to_string(ratio) << ", Warning!" << std::endl;
    }
}

// Get quadrature properties background volume cell 
std::vector<GQTriangleRPIM::CellProperties>
    PointcloudRPIM::get_quadrature_volume_cells_properties(void)
{
    return m_volume_cell_props;
}

// Get quadrature weights of background volume cell
std::vector<double> PointcloudRPIM::get_quadrature_volume_cells_weights(void)
{
    return m_volume_cell_quadr_weights;
}

// Get quadrature properties background surface cell 
std::vector<GQLineRPIM::CellProperties>
    PointcloudRPIM::get_quadrature_surface_cells_properties(void)
{
    return m_surface_cells_props;
}

// Get quadrature weights of background surface cell
std::vector<double> PointcloudRPIM::get_quadrature_surface_cells_weights(void)
{
    return m_surface_cell_quadr_weights;
}