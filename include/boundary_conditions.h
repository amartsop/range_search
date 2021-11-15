#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>

#include "geom.h"
#include "gp_utils.h"
#include "mesh2D.h"


class BoundaryConditions
{
    public:
        BoundaryConditions(const Mesh2D& mesh_raw);

        // Get mesh
        Mesh2D get_mesh(void) { return m_mesh; }

        // Get number of boundaries 
        size_t get_number_of_boundaries(void) { return m_boundaries_num; }

    private: 

        // Raw mesh
        Mesh2D m_mesh_raw;

        // Triangular mesh
        Mesh2D m_mesh;

        // Triangular mesh initial
        Mesh2D m_mesh_initial;

    private:
        // Boundary x scale factor
        double m_bound_x_scale_factor = 0.02;

        // Set boundary point indices
        void set_boundary_pts_indices(void);

        // Boundary points indices
        std::vector<size_t> m_boundary_pts_indices;

        // Number of boundary points
        size_t m_boundaries_num;

        // Free point indices (non boundaries)
        std::vector<size_t> m_free_pts_indices;

        // Number of free points 
        size_t m_free_pts_num;

        // Update boundaries
        void update_boundaries(void);

        // Generate mesh for field notes
        void generate_field_nodes_mesh(void);

        // Get the index of a specified element in a vector
        size_t get_index_of_specified_value(std::vector<size_t> v,
            size_t k);
};