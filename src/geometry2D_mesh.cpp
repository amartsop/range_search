#include "../include/geometry2D_mesh.h"

Geometry2DMesh::Geometry2DMesh(const std::string& filename2D, unsigned int boundary_density)
{
    // Read matrix from file
    m_input_mat.load(filename2D, arma::csv_ascii);

    // Generate source boundary points container
    auto source_boundary_points = generate_input_points(m_input_mat.col(0),
        m_input_mat.col(1), m_input_mat.col(2));

    // Generate adjusted points container
    auto adjusted_boundary_points = boundary_density_adjustment(source_boundary_points,
        boundary_density);
    
    // Generate mesh
    generate_mesh(adjusted_boundary_points);
}


// Generate input points 
geom::PointCloud<double> Geometry2DMesh::generate_input_points(const arma::dvec& x,
        const arma::dvec& y, const arma::dvec& z)
{
    // Get number of points
    unsigned int point_num = x.n_rows;

    // Initialize pointcloud
    geom::PointCloud<double> pc;

    for (size_t i = 0; i < point_num; i++)
    {
        // Generate point i
        geom::Point<double> point_i;

        // Point point_i;
        point_i.x = x.at(i); point_i.y = y.at(i); point_i.z = z.at(i);

        // Push back points
        pc.pts.push_back(point_i);
    }

    return pc;
}

// Gives a vector of Points which are sampled as equally-spaced segments
// taken along the linear interpolation between points in the source.
geom::PointCloud<double> Geometry2DMesh::boundary_density_adjustment(const
    geom::PointCloud<double>& source_pc, std::size_t target_count)
{
    // Initialize output
    geom::PointCloud<double> result;

    // Define source vector of points
    auto source = source_pc.pts;

    if(source.size() < 2 || target_count < 2) {
        // degenerate source vector or target_count value
        // for simplicity, this returns an empty result
        // but special cases may be handled when appropriate for the application
        return result;
    }

    // total_length is the total length along a linear interpolation
    // of the source points.
    const double total_length = linear_curve_length(source);

    // segment_length is the length between result points, taken as
    // distance traveled between these points on a linear interpolation
    // of the source points.  The actual Euclidean distance between
    // points in the result vector can vary, and is always less than
    // or equal to segment_length.
    const double segment_length = total_length / (target_count - 1);

    // start and finish are the current source segment's endpoints
    auto start = source.begin();
    auto finish = start + 1;

    // src_segment_offset is the distance along a linear interpolation
    // of the source curve from its first point to the start of the current
    // source segment.
    double src_segment_offset = 0;

    // src_segment_length is the length of a line connecting the current
    // source segment's start and finish points.
    double src_segment_length = euclidean_distance(*start, *finish);

    // The first point in the result is the same as the first point
    // in the source.
    result.pts.push_back(*start);

    for(std::size_t i=1; i<target_count-1; ++i) {
        // next_offset is the distance along a linear interpolation
        // of the source curve from its beginning to the location
        // of the i'th point in the result.
        // segment_length is multiplied by i here because iteratively
        // adding segment_length could accumulate error.
        const double next_offset = segment_length * i;

        // Check if next_offset lies inside the current source segment.
        // If not, move to the next source segment and update the
        // source segment offset and length variables.
        while(src_segment_offset + src_segment_length < next_offset) {
            src_segment_offset += src_segment_length;
            start = finish++;
            src_segment_length = euclidean_distance(*start, *finish);
        }
        // part_offset is the distance into the current source segment
        // associated with the i'th point's offset.
        const double part_offset = next_offset - src_segment_offset;

        // part_ratio is part_offset's normalized distance into the 
        // source segment. Its value is between 0 and 1 ,
        // where 0 locates the next point at "start" and1
        // locates it at "finish".  In-between values represent a
        // weighted location between these two extremes.
        const double part_ratio = part_offset / src_segment_length;

        // Use part_ratio to calculate the next point's components
        // as weighted averages of components of the current
        // source segment's points.
        result.pts.push_back({
            start->x + part_ratio * (finish->x - start->x),
            start->y + part_ratio * (finish->y - start->y),
            0.0
        });
    }

    // The first and last points of the result are exactly
    // the same as the first and last points from the input,
    // so the iterated calculation above skips calculating
    // the last point in the result, which is instead copied
    // directly from the source vector here.
    // result.push_back(source.back());
    result.pts.push_back(source.back());

    return result;
}

// Calculates the total length of the linear interpolation through a vector of Points.
double Geometry2DMesh::linear_curve_length(const std::vector<geom::Point<double>>& points)
{
    auto start = points.begin();

    if(start == points.end()) return 0;

    auto finish = start + 1;
    double sum = 0;
    while(finish != points.end()) {
        sum += euclidean_distance(*start, *finish);
        start = finish++;
    }
    return sum;
}

// Euclidean distance between two points
double Geometry2DMesh::euclidean_distance(const geom::Point<double>& a, const
        geom::Point<double>& b)
{
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double lsq = dx*dx + dy*dy;
    return std::sqrt(lsq);
}


// Generate mesh
void Geometry2DMesh::generate_mesh(const geom::PointCloud<double>& boundary_points)
{
    /*********************** Generate mesh **************************/
    // Mesh generation
    gmsh::initialize();
    gmsh::model::add("t1");

    // Points 
    std::vector<int> points(boundary_points.pts.size());

    for (size_t i = 0; i < points.size(); i++)
    {
        points.at(i) = gmsh::model::geo::addPoint(boundary_points.pts.at(i).x, 
            boundary_points.pts.at(i).y, boundary_points.pts.at(i).z);
    }

    // Lines 
    std::vector<int> lines(boundary_points.pts.size());

    // Points identifiers (current pi, next pn)
    int pi, pn;

    for (size_t i = 0; i < lines.size(); i++)
    {
        // Current point
        pi = points.at(i);

        // Next point
        pn = (i == lines.size()-1) ? points.at(0) : points.at(i+1);

        lines.at(i) = gmsh::model::geo::addLine(pi, pn);
    }

    // Add curve loop
    int curve = gmsh::model::geo::addCurveLoop(lines);

    // Add surface
    int surface = gmsh::model::geo::addPlaneSurface({curve});

    gmsh::model::geo::synchronize();

    // We can then generate a 2D mesh...
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

    // Get nodes number
    m_mesh.nodes_num = node_tags.size();

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
        LineElement line_element_i;
        line_element_i.node1_index = element_node_tags.at(line_index).at(i) - 1;
        line_element_i.node2_index = element_node_tags.at(line_index).at(i+1) - 1;

        // Get indices of boundary nodes
        m_mesh.bound_nodes_indices.push_back(line_element_i.node1_index);
        m_mesh.bound_nodes_indices.push_back(line_element_i.node2_index);

        // Push results to line element vector
        m_mesh.bound_elements.push_back(line_element_i);
    }

    // Find unique boundary indices
    std::sort(m_mesh.bound_nodes_indices.begin(), m_mesh.bound_nodes_indices.end());
    m_mesh.bound_nodes_indices.resize(std::distance(m_mesh.bound_nodes_indices.begin(), 
        std::unique(m_mesh.bound_nodes_indices.begin(), m_mesh.bound_nodes_indices.end())));  

    /******************* Triangular Elements *********************/
    // Assign triangle element tags and element num
    m_mesh.element_tags = element_tags.at(triangle_index);
    m_mesh.elements_num = element_tags.at(triangle_index).size();

    int max_index = 0;

    // Element node tags
    for (size_t i = 0; i < element_node_tags.at(triangle_index).size(); i=i+m_triangle_nodes)
    {
        Triangular2DElement element_i;
        element_i.node1_index = element_node_tags.at(triangle_index).at(i) - 1;
        element_i.node2_index = element_node_tags.at(triangle_index).at(i+1) - 1;
        element_i.node3_index = element_node_tags.at(triangle_index).at(i+2) - 1;

        if (element_i.node1_index >= max_index) { max_index = element_i.node1_index; }
        if (element_i.node2_index >= max_index) { max_index = element_i.node2_index; }
        if (element_i.node3_index >= max_index) { max_index = element_i.node3_index; }

        // Get indices of all nodes
        m_mesh.inner_nodes_indices.push_back(element_i.node1_index);
        m_mesh.inner_nodes_indices.push_back(element_i.node2_index);
        m_mesh.inner_nodes_indices.push_back(element_i.node3_index);

        // Push results to elements vector
        m_mesh.elements.push_back(element_i);
    }

    // Find unique inner indices
    std::sort(m_mesh.inner_nodes_indices.begin(), m_mesh.inner_nodes_indices.end());
    m_mesh.inner_nodes_indices.resize(std::distance(m_mesh.inner_nodes_indices.begin(), 
        std::unique(m_mesh.inner_nodes_indices.begin(), m_mesh.inner_nodes_indices.end())));  

    // Remove boundary indices
    for(size_t i = 0; i < m_mesh.bound_nodes_indices.size(); i++)
    {
        m_mesh.inner_nodes_indices.erase(std::remove(
            m_mesh.inner_nodes_indices.begin(), m_mesh.inner_nodes_indices.end(),
            m_mesh.bound_nodes_indices.at(i)), m_mesh.inner_nodes_indices.end());
    }
}

// Get boundary points
geom::PointCloud<double> Geometry2DMesh::get_boundary_points(const
    Geometry2DMesh::Triangular2DMesh & mesh)
{
    // Initialize pointcloud
    geom::PointCloud<double> pc;

    for (size_t i = 0; i < mesh.bound_nodes_indices.size(); i++)
    {
        // Get index
        size_t idx = mesh.bound_nodes_indices.at(i);

        // std::cout << idx << std::endl;

        // Initialize point and assign to stuct
        geom::Point<double> point_i;
        point_i.x = mesh.node_coords.pts.at(idx).x;
        point_i.y = mesh.node_coords.pts.at(idx).y;
        point_i.z = mesh.node_coords.pts.at(idx).z;

        // Push to points
        pc.pts.push_back(point_i);
    }

    return  pc;
}

// Get inner points
geom::PointCloud<double> Geometry2DMesh::get_inner_points(const
    Geometry2DMesh::Triangular2DMesh & mesh)
{
    // Initialize pointcloud
    geom::PointCloud<double> pc;

    for (size_t i = 0; i < mesh.inner_nodes_indices.size(); i++)
    {
        // Get index
        size_t idx = mesh.inner_nodes_indices.at(i);

        // Initialize point and assign to stuct
        geom::Point<double> point_i;
        point_i.x = mesh.node_coords.pts.at(idx).x;
        point_i.y = mesh.node_coords.pts.at(idx).y;
        point_i.z = mesh.node_coords.pts.at(idx).z;


        // Push to points
        pc.pts.push_back(point_i);
    }

    return  pc;
}