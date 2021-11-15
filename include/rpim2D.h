#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>

#include "geom.h"
#include "gp_utils.h"
#include "./boundary_conditions.h"
#include "./pointcloud_rpim.h"

#include "./mesh2D.h"

#include "./gq_triangle_rpim.h"
#include "./gq_line_rpim.h"
#include "./support_domain.h"
#include "./geometry_model.h"
#include "./shape_function.h"
#include "./strain.h"
#include "./point_loads.h"

#include <boost/tuple/tuple.hpp>
#include "gnuplot-iostream.h"
#include <bits/stdc++.h>

class RPIM2D 
{
    public:
        RPIM2D() {}

        // Initialize 
        void initialize(const Mesh2D& mesh_raw, double thickness,
            const geom::RPIMParameters& params, bool animate=false);

        // Update field nodes
        void update(const arma::dvec& q_bar);

        // Get number of dofs
        size_t get_dofs_num(void) { return m_dofs_num; }

        // Get initial field nodes mesh
        Mesh2D get_initial_field_nodes_mesh(void) { return m_field_nodes_mesh_initial; }

    public:

        // Active external force vector getter
        arma::dvec get_external_force_vector_a(void) const { return m_fex_a; }

        // Constrained external force vector getter
        arma::dvec get_external_force_vector_c(void) const { return m_fex_c; }

        // Get kcc matrix
        arma::dmat get_kcc_matrix(void) const { return m_kfcc; }

        // Get kca matrix
        arma::dmat get_kca_matrix(void) const { return m_kfca; }

        // Get kaa matrix
        arma::dmat get_kaa_matrix(void) const { return m_kfaa; }

    public:

        // Get full state vector
        arma::dvec get_full_state_vector(const arma::dvec& qa);

        // Get boundary state vector
        arma::dvec get_boundary_state_vector(const arma::dvec& q_bar);

        // Get deformed state of pointcloud of interest
        geom::PointCloud<double> get_deformed_state(const
            geom::PointCloud<double>& inter_pc, const arma::dvec& q_bar) const;

        // Get deformed mesh
        Mesh2D get_deformed_mesh(const arma::dvec& q_bar) const;
        
    private:

        // Mesh
        Mesh2D m_field_nodes_mesh, m_field_nodes_mesh_initial;
        
        // Number of field nodes
        size_t m_field_nodes_num;

        // Dofs per node
        int m_dofs_per_node = 2;

        // Total number of dofs
        size_t m_dofs_num;

        // Number of boundaries
        size_t m_boundaries_num;

        // Thickness (m)
        double m_thickness;

    private:
        // External forces
        arma::dvec m_fex_c, m_fex_a;

        // Stiffness matrices
        arma::dmat m_kfcc, m_kfca, m_kfaa;

        // Point loads
        PointLoads m_point_load;

    private:        

        // Pointcloud handler 
        PointcloudRPIM m_pc_rpim;

        // Update field nodes mesh and cloud    
        void update_field_nodes_mesh_and_cloud(const arma::dvec& q_bar);

        // Global cloud
        geom::PointCloud<double> m_cloud, m_cloud_initial;

        // Support domain handle
        SupportDomain m_sup_domain;
    
        // Support domain radius 
        geom::RPIMParameters m_search_params;

        // Support domain animation flag
        bool m_animation_flag;
    
    private:

        // Quadrature sampling for surface elements
        const short int m_volume_quadr_interp = GQTriangleRPIM::QUADRATIC;

        // Quadrature properties of background volume cell
        std::vector<GQTriangleRPIM::CellProperties> m_volume_cell_props;

        // Quadrature weights of background volume cell
        std::vector<double> m_volume_cell_quadr_weights;

    private: 

        // Quadrature sampling for surface elements
        const short int m_surface_quadr_interp = 3; 

        // Quadrature properties of background surface cell
        std::vector<GQLineRPIM::CellProperties> m_surface_cell_props;

        // Quadrature weights of background surface cell
        std::vector<double> m_surface_cell_quadr_weights;
};