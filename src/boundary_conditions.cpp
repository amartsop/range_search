#include "../include/boundary_conditions.h"

BoundaryConditions::BoundaryConditions(const Mesh2D& mesh_raw)
{
    // Get mesh (reference mesh)
    m_mesh_raw = mesh_raw;

    // Set boundary points indices
    set_boundary_pts_indices();

    // Generate field nodes mesh
    generate_field_nodes_mesh();

    // Update boundary indices
    update_boundaries();
}

// Set boundary point indices
void BoundaryConditions::set_boundary_pts_indices(void)
{
    /*Set the boundaries as the points that are less than x_thresh of the x
    range */
    // Get field nodes vector
    auto field_nodes_vec = geom::conv_pc_to_pc_vec<double>(
        m_mesh_raw.node_coords);

    // Calculate range in y direction
    auto range = gp_utils::find_vector_range<double>(field_nodes_vec.x);

    // Define threshold in y direction
    double x_thresh = range.min + m_bound_x_scale_factor * std::abs(range.val);

    // Loop through field nodes and assign boundary indices and pointsj
    for (size_t i = 0; i < m_mesh_raw.node_coords.pts.size(); i++)
    {
        if (m_mesh_raw.node_coords.pts.at(i).x < x_thresh)
        {
            m_boundary_pts_indices.push_back(i);
        }
    }
    
    // Assign number of boundary points
    m_boundaries_num = m_boundary_pts_indices.size();

    // Create free points indices
    for(size_t i = 0; i < m_mesh_raw.node_coords.pts.size(); i++)
    {
        m_free_pts_indices.push_back(i);
    }

    for (size_t i = 0; i < m_boundary_pts_indices.size(); i++)
    {
        m_free_pts_indices.erase(m_free_pts_indices.begin() +
            m_boundary_pts_indices.at(i)-i);
    }

    // Assign number of free points
    m_free_pts_num = m_free_pts_indices.size();
}


// Generate field nodes mesh
void BoundaryConditions::generate_field_nodes_mesh(void)
{
    // Indices of field nodes mesh
    m_mesh.node_indices = m_mesh_raw.node_indices;

    // Desired indices format
    std::vector<size_t> state_format;
    state_format.insert(state_format.end(), m_boundary_pts_indices.begin(), 
        m_boundary_pts_indices.end());
    state_format.insert(state_format.end(), m_free_pts_indices.begin(), 
        m_free_pts_indices.end());

    // Nodes coordinates for field nodes mesh
    for (size_t i = 0; i < state_format.size(); i++)
    {
        auto pt_i = m_mesh_raw.node_coords.pts.at(state_format.at(i));
        m_mesh.node_coords.pts.push_back(pt_i);
    }

    // Elements for field nodes mesh
    for (size_t i = 0; i < m_mesh_raw.volume_elements.size(); i++)
    {
        // Element i for triangle mesh
        auto element_i_triangle_mesh = m_mesh_raw.volume_elements.at(i);

        // Initialize element i for field mesh
        geom::Triangular2DElement element_i_field_mesh;
        
        // Find index 
        element_i_field_mesh.node1_index =
            get_index_of_specified_value(state_format, 
            element_i_triangle_mesh.node1_index);

        element_i_field_mesh.node2_index =
            get_index_of_specified_value(state_format, 
            element_i_triangle_mesh.node2_index);
        
        element_i_field_mesh.node3_index =
            get_index_of_specified_value(state_format, 
            element_i_triangle_mesh.node3_index);

        // Push to field nodes mesh
        m_mesh.volume_elements.push_back(element_i_field_mesh);
    }

    // Boundary elements for field nodes mesh
    for (size_t i = 0; i < m_mesh_raw.bound_elements.size(); i++)
    {
        // Boundaery element i for triangle mesh
        auto b_element_i_triangle_mesh = m_mesh_raw.bound_elements.at(i);

        // Initialize boundary element i for field mesh
        geom::LineElement b_element_i_field_mesh;

        // Assign indices
        b_element_i_field_mesh.node1_index = get_index_of_specified_value(
            state_format, b_element_i_triangle_mesh.node1_index);

        b_element_i_field_mesh.node2_index = get_index_of_specified_value(
            state_format, b_element_i_triangle_mesh.node2_index);

        // Push to field nodes mesh
        m_mesh.bound_elements.push_back(b_element_i_field_mesh);
    }

    // Initial field nodes mesh
    m_mesh_initial = m_mesh;
}

// Get the index of a specified element in a vector
size_t BoundaryConditions::get_index_of_specified_value(
    std::vector<size_t> v, size_t k)
{
    auto it = std::find(v.begin(), v.end(), k);
 
    int idx=-1;

    // If element was found
    if (it != v.end()) { idx = it - v.begin(); }

    return idx;
}


// Update boundaries
void BoundaryConditions::update_boundaries(void)
{
    // Update boundary points
    for (size_t i = 0; i < m_boundary_pts_indices.size(); i++)
    {
        m_boundary_pts_indices.at(i) = i;
    }

    // Update free points
    for (size_t i = 0; i < m_free_pts_indices.size(); i++)
    {
        m_free_pts_indices.at(i) = i + m_boundary_pts_indices.size();
    }
}