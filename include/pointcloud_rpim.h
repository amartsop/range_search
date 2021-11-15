#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>

#include "geom.h"
#include "gp_utils.h"
#include "./mesh2D.h"
#include "./gq_triangle_rpim.h"
#include "./gq_line_rpim.h"


class PointcloudRPIM
{
    public:
        PointcloudRPIM() {};
        
        // Initialize pointcloud
        void initialize(const Mesh2D& field_nodes_mesh,
            size_t volume_quadr_interp, size_t surface_quadr_interp);

        // Get global cloud
        geom::PointCloud<double> get_cloud(void) { return m_cloud; }

        // Get number of volume quadrature points
        size_t get_number_of_volume_quadrature_points(void) { return m_volume_gq_pts_num; }

    public:
        // Get number of volume cells
        size_t get_volume_cells_number(void) { return m_volume_cell_props.size(); }

        // Get quadrature properties background volume cell 
        std::vector<GQTriangleRPIM::CellProperties>
            get_quadrature_volume_cells_properties(void);

        // Get quadrature weights of background volume cell
        std::vector<double> get_quadrature_volume_cells_weights(void);

    public:

        // Get number of surface cells
        size_t get_surface_cells_number(void) { return m_surface_cells_props.size(); }

        // Get quadrature properties background surface cell 
        std::vector<GQLineRPIM::CellProperties>
            get_quadrature_surface_cells_properties(void);

        // Get quadrature weights of background surface cell
        std::vector<double> get_quadrature_surface_cells_weights(void);
        
    private: 

        // Mesh
        Mesh2D m_field_nodes_mesh;

        // Number of field nodes
        size_t m_field_nodes_num;

        // Generate data points structure
        void generate_data_points(void);

        // Quadrature point check
        void quadrature_points_check(void);

        // Global cloud
        geom::PointCloud<double> m_cloud;

    private:
        // Generate volume quadrature points
        void generate_volume_quadrature_points(void);

        // Volume quadrature points
        geom::PointCloud<double> m_volume_gq_pts;

        // Number of volume quadrature points
        size_t m_volume_gq_pts_num;

        // Volume quadrature handle
        GQTriangleRPIM m_volume_gq;

        // Quadrature sampling for volume element
        size_t m_volume_quadr_interp;

        // Background volume cell properties
        std::vector<GQTriangleRPIM::CellProperties> m_volume_cell_props;

        // Weights of background volume cell
        std::vector<double> m_volume_cell_quadr_weights;

    private:

        // Generate surface quadrature points
        void generate_surface_quadrature_points(void);

        // Surface quadrature points
        geom::PointCloud<double> m_surface_gq_pts;

        // Number of surface quadrature points
        size_t m_surface_gq_pts_num;

        // Surface quadrature handle
        GQLineRPIM m_surface_gq;

        // Quadrature sampling for surface elements
        size_t m_surface_quadr_interp; 
        
        // Background surface cell properties
        std::vector<GQLineRPIM::CellProperties> m_surface_cells_props;

        // Weights of background surface cell
        std::vector<double> m_surface_cell_quadr_weights;
};