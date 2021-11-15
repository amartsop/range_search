#pragma once 

#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <armadillo>
#include <algorithm>

#include <boost/tuple/tuple.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_2.h>

#include "geom.h"
#include "gp_utils.h"
#include "gnuplot-iostream.h"

#include "kd_trees.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_set_2<K>::Vertex_handle Vertex_handle;


class SupportDomain
{

public:
    SupportDomain() {};

    struct SupportDomainPoint
    {
        // Point idx (that's a gaussian pointgTgt)
        size_t point_idx;

        // Coordinates of index point
        arma::dvec point_coords;

        // Supporting nodes indices
        std::vector<size_t> support_indices;

        // Coordinates of support_indices
        std::vector<arma::dvec> support_coords;
    };
    
    // Generate support domain
    std::vector<SupportDomainPoint> generate(
        const geom::PointCloud<double>& field_nodes,
        const geom::PointCloud<double>& data_pts,
        const geom::RPIMParameters& rpim_params, bool animate=false);

private:

    // Pointcloud to GpointCloud
    std::vector<K::Point_2> conv_pc_to_gpc(const geom::PointCloud<double>& pc);

    // Rectange structure
    geom::PointCloudVec<double> rectangle(double x_center, double y_center,
        double width, double height);

    // Get equal index
    size_t get_equal_idx(const K::Point_2& inter_point, const
        std::vector<K::Point_2>& data_pts, double tol = 1.0e-10);

private:

    // Gnuplot handle 
    Gnuplot m_gp;

    // Animate support domain
    void animate_support_domain(const geom::PointCloud<double>& cloud,
        const geom::PointCloud<double>& field_nodes,
        const std::vector<SupportDomainPoint>& sup_dom_pts, double rect_width,
        double rect_height);
};