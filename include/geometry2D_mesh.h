#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>
#include <algorithm>
#include <iterator>

#include "geom.h"

#include <set>
#include <gmsh.h>

#include "gnuplot-iostream.h"
#include"gp_utils.h"


class Geometry2DMesh
{

public:
    Geometry2DMesh(const std::string& filename2D, unsigned int boundary_density);

// Triangular 2D element struct
struct Triangular2DElement
{
    int node1_index;
    int node2_index;
    int node3_index;
};

// Boundary line element struct
struct LineElement
{
    int node1_index;
    int node2_index;
};

// Triangular 2D mesh struct
struct Triangular2DMesh {

    // Number of nodes
    int nodes_num;

    // Mesh nodes indices
    std::vector<std::size_t> node_indices;

    // Coordinates of mesh nodes
    geom::PointCloud<double> node_coords;

    // Number of elements
    int elements_num;

    // Elements tags
    std::vector<std::size_t> element_tags;

    // Elements
    std::vector<Triangular2DElement> elements;

    // Boundary elements
    std::vector<LineElement> bound_elements;

    // Inner points indices
    std::vector<std::size_t> inner_nodes_indices;

    // Boundary points indices
    std::vector<std::size_t> bound_nodes_indices;
};

// Get boundary points
static geom::PointCloud<double> get_boundary_points(const
    Geometry2DMesh::Triangular2DMesh & mesh);

// Get inner points
static geom::PointCloud<double> get_inner_points(const
    Geometry2DMesh::Triangular2DMesh & mesh);

// Get triangular 2D mesh
Triangular2DMesh get_triangular2D_mesh(void) { return m_mesh; }

private:

    // Input matrix
    arma::dmat m_input_mat;
    
    // Boundary and inner points
    geom::PointCloud<double> m_bound;
    geom::PointCloud<double> m_inner;

    // Initialize mesh 
    Triangular2DMesh m_mesh;

    // Generate input points 
    geom::PointCloud<double> generate_input_points(const arma::dvec& x,
        const arma::dvec& y, const arma::dvec& z);

    // Adjust point density
    geom::PointCloud<double> boundary_density_adjustment(const
        geom::PointCloud<double>& source, std::size_t target_count);

    // Calculates the total length of the linear interpolation through a vector of Points.
    double linear_curve_length(const std::vector<geom::Point<double>>& points);

    // Euclidean distance between two points
    double euclidean_distance(const geom::Point<double>& a, const
        geom::Point<double>& b);

    // Generate mesh
    void generate_mesh(const geom::PointCloud<double>& boundary_points);

    // 2D Triangle id and number of nodes per triangle (GMSH constants)
    const int m_triangle2D_id = 2;
    const int m_triangle_nodes = 3;

    // Boundary line id and number of nodes per line (GMSH consants)
    const int m_line_id = 1;
    const int m_line_nodes = 2;
};