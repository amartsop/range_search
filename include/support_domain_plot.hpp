#pragma once 

#include <iostream>
#include <vector>
#include "geom.h"
#include "gnuplot-iostream.h"
#include "geometry2D_mesh.h"

class SupportDomainPlot
{
public:
    SupportDomainPlot(const geom::SupportDomain& sup_dom,
        const Geometry2DMesh::Triangular2DMesh& mesh);

    // Update plot
    void update(const geom::Point<double>& query_pt, const std::vector<int>& 
        sup_dom_indices);

private:

    // Support domain properties
    std::vector<double> m_circle_x{0.0};
    std::vector<double> m_circle_y{0.0};
    std::vector<double> m_circle_r{0.0};

    // Mesh member variable
    Geometry2DMesh::Triangular2DMesh m_mesh;

    // Inner and boundary points
    geom::PointCloud<double> m_bound_pts, m_inner_pts;

    // Interest point
    std::vector<double> m_quer_x{0.0};
    std::vector<double> m_quer_y{0.0};

    // Gnuplot handle
    Gnuplot m_gp;
};

SupportDomainPlot::SupportDomainPlot(const geom::SupportDomain& sup_dom,
        const Geometry2DMesh::Triangular2DMesh& mesh)
{
    // Initialize support domain propertiesj
    m_circle_x.at(0) = sup_dom.xc; m_circle_y.at(0) = sup_dom.xc;
    m_circle_r.at(0) = sup_dom.r;

    // Store pointcloud
    m_mesh = mesh;
    
    // Initialize plos
    m_gp << "set xrange [-12.0:12.0]\nset yrange [-12.0:12.0]\n";
    // m_gp << "set xrange [-1.5:1.5]\nset yrange [-1.5:1.5]\n";
    m_gp << "set style fill transparent solid 0.1 \n";

    // Get inner and boundary points
    m_bound_pts = Geometry2DMesh::get_boundary_points(m_mesh);
    m_inner_pts = Geometry2DMesh::get_inner_points(m_mesh);
}

// Update plot
void SupportDomainPlot::update(const geom::Point<double>& query_pt,
    const std::vector<int>& sup_dom_indices)
{
    // Update support domain
    m_circle_x.at(0) = query_pt.x;
    m_circle_y.at(0) = query_pt.y;

    // Update query point
    m_quer_x.at(0) = query_pt.x;
    m_quer_y.at(0) = query_pt.y;

    std::vector<double> sup_dom_x, sup_dom_y;
    
    // Get support domain points
    for (size_t i = 0; i < sup_dom_indices.size(); i++)
    {
        // Get index
        size_t idx = sup_dom_indices.at(i);

        // Get point
        sup_dom_x.push_back(m_mesh.node_coords.pts.at(idx).x);
        sup_dom_y.push_back(m_mesh.node_coords.pts.at(idx).y);
    }

    // Convert pointcloud to pointcloud vector
    geom::PointCloudVec<double> inner_pts_vec = geom::conv_pc_to_pc_vec<double>(m_inner_pts);
    geom::PointCloudVec<double> bound_pts_vec = geom::conv_pc_to_pc_vec<double>(m_bound_pts);


    // Animate
    std::string plot_str = "plot '-' with points pt 7 ps 0.5 lc rgb 'red', "
        "'-' with points pt 2 ps 1.0 lc rgb 'blue', "
        "'-' with points pt 1 ps 1.0 lc rgb 'black', "
        "'-' with lines, '-' with circles\n";
    
    m_gp << plot_str;
    
    // Send data
    m_gp.send1d(boost::make_tuple(inner_pts_vec.x, inner_pts_vec.y));
    m_gp.send1d(boost::make_tuple(m_quer_x, m_quer_y));
    m_gp.send1d(boost::make_tuple(sup_dom_x, sup_dom_y));
    m_gp.send1d(boost::make_tuple(bound_pts_vec.x, bound_pts_vec.y));
    m_gp.send1d(boost::make_tuple(m_circle_x, m_circle_y, m_circle_r));
}