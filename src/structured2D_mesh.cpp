#include "../include/structured2D_mesh.h"

Structured2DMesh::Structured2DMesh(const geom::Point<double>& pt1,
    const geom::Point<double>& pt2, int pts_x, int pts_y)
{
    // Get points in x and y       
    m_points_x = pts_x; m_points_y = pts_y;

    // Initialize mesh
    gmsh::initialize();

    gmsh::model::add("t1");

    // Tolerance
    double lc = 1e-2;

    // Generate points
    int p1 = gmsh::model::geo::addPoint(pt1.x, pt1.y, pt1.z, lc);
    int p2 = gmsh::model::geo::addPoint(pt1.x, pt2.y, pt1.z, lc);
    int p3 = gmsh::model::geo::addPoint(pt2.x, pt2.y, pt1.z, lc);
    int p4 = gmsh::model::geo::addPoint(pt2.x, pt1.y, pt1.z, lc);

    // Add lines
    int l1 = gmsh::model::geo::addLine(p1, p2);
    int l2 = gmsh::model::geo::addLine(p2, p3);
    int l3 = gmsh::model::geo::addLine(p3, p4);
    int l4 = gmsh::model::geo::addLine(p4, p1);

    // Add curve loops
    std::vector<int> lines = {l1, l2, l3, l4};
    int curve = gmsh::model::geo::addCurveLoop(lines);

    // Add surface
    int surface = gmsh::model::geo::addPlaneSurface({curve});

    // Set spacing
    gmsh::model::geo::mesh::setTransfiniteCurve(l1, pts_y);
    gmsh::model::geo::mesh::setTransfiniteCurve(l2, pts_x);
    gmsh::model::geo::mesh::setTransfiniteCurve(l3, pts_y);
    gmsh::model::geo::mesh::setTransfiniteCurve(l4, pts_x);

    gmsh::model::geo::mesh::setTransfiniteSurface(surface, "Left", lines);

    // Smoothing
    gmsh::option::setNumber("Mesh.Smoothing", 100);

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(2);

    /*********************** Get mesh properties **************************/

    /***************** Nodes *********************/
    // Initialize node vectors
    std::vector<size_t> node_tags;
    std::vector<double> node_coords;
    std::vector<double> parametric_coords;

    // Get node properties
    gmsh::model::mesh::getNodes(node_tags, node_coords, parametric_coords, -1,
        -1, false, false);

    // Get node indices
    for (size_t i = 0; i < node_tags.size(); i++)
    {
        m_mesh.node_indices.push_back(node_tags.at(i) - 1);
    }

    // Get node coordinates
    for(size_t i = 0; i < node_coords.size(); i=i+m_triangle_nodes)
    {
        // Initialize points
        geom::Point<double> point_i;
        point_i.x = node_coords.at(i);
        point_i.y = node_coords.at(i+1);
        point_i.z = node_coords.at(i+2);

        // Push results to node_coords vector
        m_mesh.node_coords.pts.push_back(point_i);
    }

    /***************** Elements *********************/

    // Get element properties
    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t> > element_node_tags;

    gmsh::model::mesh::getElements(element_types, element_tags, element_node_tags);

    // Get the index of 2D triangles and boundary lines
    int triangle_index = 0;
    int line_index = 0;
    for (size_t i = 0; i < element_types.size(); i++)
    {
        if (element_types.at(i) == m_triangle2D_id){ triangle_index = i; }
        if (element_types.at(i) == m_line_id){ line_index = i; }
    }

    /******************* Line Elements at Boundary *********************/
    // Boundary element node tags
    for (size_t i = 0; i < element_node_tags.at(line_index).size(); i=i+m_line_nodes)
    {
        // Line elements
        geom::LineElement line_element_i;
        line_element_i.node1_index = element_node_tags.at(line_index).at(i) - 1;
        line_element_i.node2_index = element_node_tags.at(line_index).at(i+1) - 1;

        // Push results to line element vector
        m_mesh.bound_elements.push_back(line_element_i);
    }


    /******************* Triangular Elements *********************/
    int max_index = 0;
    // Element node tags
    for (size_t i = 0; i < element_node_tags.at(triangle_index).size(); i=i+m_triangle_nodes)
    {
        geom::Triangular2DElement element_i;
        element_i.node1_index = element_node_tags.at(triangle_index).at(i) - 1;
        element_i.node2_index = element_node_tags.at(triangle_index).at(i+1) - 1;
        element_i.node3_index = element_node_tags.at(triangle_index).at(i+2) - 1;

        if (element_i.node1_index >= max_index) { max_index = element_i.node1_index; }
        if (element_i.node2_index >= max_index) { max_index = element_i.node2_index; }
        if (element_i.node3_index >= max_index) { max_index = element_i.node3_index; }

        // Push results to elements vector
        m_mesh.volume_elements.push_back(element_i);
    }

    // Calculate nodal spacing
    m_dcx = std::abs(pt2.x - pt1.x) / (pts_x-1);
    m_dcy = std::abs(pt2.y - pt1.y) / (pts_y-1);

    m_dc = std::sqrt(pow(m_dcx, 2.0) + pow(m_dcy, 2.0));
}
