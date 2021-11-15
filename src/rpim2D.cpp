#include "../include/rpim2D.h"

void RPIM2D::initialize(const Mesh2D& mesh_raw, double thickness,
    const geom::RPIMParameters& params, bool animate)
{
    // Set search params 
    m_search_params = params;

    // Set animation flag
    m_animation_flag = animate;

    // Initialize boundary conditions
    BoundaryConditions bc(mesh_raw);

    // Geometry thickness
    m_thickness = thickness;

    // Get field nodes mesh mesh
    m_field_nodes_mesh = bc.get_mesh();
    m_field_nodes_mesh_initial = m_field_nodes_mesh;
    
    // Get number of boundaries
    m_boundaries_num = bc.get_number_of_boundaries();

    // Set field nodes num
    m_field_nodes_num = m_field_nodes_mesh.node_coords.pts.size();

    // Define number of dofs
    m_dofs_num = m_dofs_per_node * m_field_nodes_num;

    // Initialize point loads
    m_point_load.initialize(m_field_nodes_mesh_initial);

    // Generate pointcloud for rpim 
    m_pc_rpim.initialize(m_field_nodes_mesh, m_volume_quadr_interp,
        m_surface_quadr_interp);

    // Get cloud
    m_cloud = m_pc_rpim.get_cloud();
    m_cloud_initial = m_cloud;

    // Set quadrature properties of background volume cell
    m_volume_cell_props = m_pc_rpim.get_quadrature_volume_cells_properties();

    // Set quadrature weights of background volume cell
    m_volume_cell_quadr_weights = m_pc_rpim.get_quadrature_volume_cells_weights();

    // Set quadrature properties of background surface cell
    m_surface_cell_props = m_pc_rpim.get_quadrature_surface_cells_properties();

    // Set quadrature weights of background volume cell
    m_surface_cell_quadr_weights = m_pc_rpim.get_quadrature_surface_cells_weights();
}

// Update field nodes
void RPIM2D::update(const arma::dvec& q_bar)
{   
    // Update field nodes mesh and cloud
    update_field_nodes_mesh_and_cloud(q_bar);

    // Generate support domain structure
    auto sup_domain_pts = m_sup_domain.generate(m_field_nodes_mesh.node_coords,
        m_cloud, m_search_params, m_animation_flag);
}    

// Update field nodes mesh
void RPIM2D::update_field_nodes_mesh_and_cloud(const arma::dvec& q_bar)
{
    // PARALLELISE
    for (size_t i = 0; i < m_field_nodes_num; i++)
    {
        // Update field nodes sturcture 
        m_field_nodes_mesh.node_coords.pts.at(i).x = 
            m_field_nodes_mesh_initial.node_coords.pts.at(i).x + q_bar(2*i);
        
        m_field_nodes_mesh.node_coords.pts.at(i).y = 
            m_field_nodes_mesh_initial.node_coords.pts.at(i).y + q_bar(2*i+1);

        // Update cloud
        m_cloud.pts.at(i).x = m_cloud_initial.pts.at(i).x + q_bar(2*i);

        m_cloud.pts.at(i).y = m_cloud_initial.pts.at(i).y + q_bar(2*i+1);
    }
}

// Get boundary state vector
arma::dvec RPIM2D::get_boundary_state_vector(const arma::dvec& q_bar)
{
    return q_bar.rows(0, 2*m_boundaries_num-1);
}

// Get deformed state of pointcloud of interest
geom::PointCloud<double> RPIM2D::get_deformed_state(const
    geom::PointCloud<double>& inter_pc, const arma::dvec& q_bar) const
{

    // Initilize final positions of interest point
    geom::PointCloud<double> final_inter_pc;  

    // Initiial point cloud
    geom::PointCloud<double> field_nodes = m_field_nodes_mesh_initial.node_coords;

    /****************** Construct data container *************************/
    // Global points number
    size_t data_pts_num = m_field_nodes_mesh_initial.node_coords.pts.size() +
        inter_pc.pts.size();

    // Reserve values and insert data points
    geom::PointCloud<double> cloud;
    cloud.pts.reserve(data_pts_num);

    // Append field nodes
    cloud.pts.insert(cloud.pts.end(), field_nodes.pts.begin(),
        field_nodes.pts.end());

    // Append interest points
    cloud.pts.insert(cloud.pts.end(), inter_pc.pts.begin(), inter_pc.pts.end());

    // Generate support domain structure
    SupportDomain sup_domain;

    auto sup_domain_pts = sup_domain.generate(field_nodes, cloud,
        m_search_params, false);

    // Initialize geometry model
    GeometryModel geom_model(sup_domain_pts, m_search_params);

    // Loop through interest points
    for (size_t i = 0; i < inter_pc.pts.size(); i++)
    {
        // Update geometric model
        geom_model.update(i, q_bar);

        // Get deformation
        auto deformation_arma = geom_model.get_deformation();

        // Calculate deformed point
        geom::Point<double> final_point;
        final_point.x = inter_pc.pts.at(i).x + deformation_arma(0);
        final_point.y = inter_pc.pts.at(i).y + deformation_arma(1);
        final_point.z = inter_pc.pts.at(i).z + 0.0;

        // Push deformed point to pc 
        final_inter_pc.pts.push_back(final_point);
    }
    
    return final_inter_pc;
}

// Get deformed mesh
Mesh2D RPIM2D::get_deformed_mesh(const arma::dvec& q_bar) const
{
    // Initialize deformed mesh
    Mesh2D deformed_mesh = m_field_nodes_mesh_initial;

    for (size_t i = 0; i < deformed_mesh.node_coords.pts.size(); i++)
    {
        deformed_mesh.node_coords.pts.at(i).x += q_bar(2*i);
        deformed_mesh.node_coords.pts.at(i).y += q_bar(2*i+1);
        deformed_mesh.node_coords.pts.at(i).z += 0.0;
    }
    return deformed_mesh;
}