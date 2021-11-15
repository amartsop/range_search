#include "../include/support_domain.h"

// Generate support domain
std::vector<SupportDomain::SupportDomainPoint> SupportDomain::generate(
    const geom::PointCloud<double>& field_nodes,
    const geom::PointCloud<double>& data_pts,
    const geom::RPIMParameters& rpim_params, bool animate)
{
    // Field nodes threshold index; any index bigger this (for the data_pts) 
    // is not a field node (aka it is a quadrature point)
    size_t field_nodes_thresh_idx = field_nodes.pts.size() - 1;

    // Initialize the vector of support domain points (that's the quadrature points)
    std::vector<SupportDomainPoint> sup_dom_pts(data_pts.pts.size() -
        field_nodes.pts.size());

    // Rectangle width and height
    double width = rpim_params.as * rpim_params.dc_x;
    double height = rpim_params.as * rpim_params.dc_y;
    
    // Define search radius
    double search_radius = std::sqrt(pow(width/2.0, 2.0) + pow(height/2.0, 2.0));

    // Initialize kd trees
    KDTrees kd_trees(field_nodes, data_pts);

    // Loop through interest points
    for (size_t query_idx = 0; query_idx < field_nodes.pts.size(); query_idx++)
    {
        // Query point
        double query_pt[3] = {field_nodes.pts.at(query_idx).x,
            field_nodes.pts.at(query_idx).y,
            field_nodes.pts.at(query_idx).z};

        std::cout << query_pt[0] << std::endl;
        std::cout << query_pt[1] << std::endl;
        std::cout << query_pt[2] << std::endl;
        std::cout << "================" << std::endl;

        // Range search
        std::vector<int> indices = kd_trees.radius_search(query_idx, search_radius);

        // Loop through found indices
        for (size_t i = 0; i < indices.size(); i++)
        {
            // Get the idx 
            auto idx = indices.at(i);

            // If the idx is for a quadrature point then add the idx of the 
            // field node to the support domain of the quadrature point
            if(idx > field_nodes_thresh_idx)
            {
                // Define goal index 
                // size_t goal_idx = idx-field_nodes.pts.size();
                size_t goal_idx = idx-field_nodes.pts.size();

                // Get the idx position of the field node
                sup_dom_pts.at(goal_idx).point_idx = idx;

                // Get the coordinates of the gaussian point 
                geom::Point<double> quadr_pt = data_pts.pts.at(idx);
                sup_dom_pts.at(goal_idx).point_coords = 
                    {quadr_pt.x, quadr_pt.y};

                // Push back the index of the field node
                sup_dom_pts.at(goal_idx).support_indices.push_back(i);

                // Push back the coordinates of the field node
                sup_dom_pts.at(goal_idx).support_coords.push_back({query_pt[0], query_pt[1]});
            }
        }
    }

    // Animate support domain
    if(animate)
    {
        animate_support_domain(data_pts, field_nodes, sup_dom_pts, width, height);
    }

    return sup_dom_pts;
}


void SupportDomain::animate_support_domain(const geom::PointCloud<double>& cloud,
    const geom::PointCloud<double>& field_nodes,
    const std::vector<SupportDomainPoint>& sup_dom_pts, double rect_width,
    double rect_height)
{
    // Initialize query point
    std::vector<double> quer_x {0.0};
    std::vector<double> quer_y {0.0};

    // Convert to vector sturcture
    geom::PointCloudVec<double> field_nodes_vec =
        geom::conv_pc_to_pc_vec<double>(field_nodes);

    geom::PointCloudVec<double> cloud_vec =
        geom::conv_pc_to_pc_vec<double>(cloud);

    // Calculate plot range
    m_gp << "set xrange " + gp_utils::plot_range(field_nodes_vec.x) + "\n";
    m_gp << "set yrange " + gp_utils::plot_range(field_nodes_vec.y) + "\n";
    m_gp << "set size ratio -1\n";

    // Loop throught quadrature points
    for (size_t i = 0; i < sup_dom_pts.size(); i++)
    {
        // Get index of quadrature point i 
        size_t quadr_pt_idx = sup_dom_pts.at(i).point_idx;

        // Get quadrature point i 
        geom::Point<double> quadr_pt = cloud.pts.at(quadr_pt_idx);

        // Query point
        quer_x.at(0) = quadr_pt.x;
        quer_y.at(0) = quadr_pt.y;

        // Rectangle
        auto rect_pc = rectangle(quadr_pt.x, quadr_pt.y, rect_width, rect_height);

        // Support domain structure
        std::vector<double> sup_dom_x, sup_dom_y;

        for (size_t k = 0; k < sup_dom_pts.at(i).support_indices.size(); k++)
        {
            // Support point (field node)
           size_t sup_idx = sup_dom_pts.at(i).support_indices.at(k);
           geom::Point<double> sup_point = cloud.pts.at(sup_idx);

            // Push support point coordinates to support vector
            sup_dom_x.push_back(sup_point.x);
            sup_dom_y.push_back(sup_point.y);
            
        }

        /************************* Animate ******************************/

        // Set linestyle 1
        // m_gp << "set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2\n";

        std::string plot_str = "plot '-' with points pt 0 ps 0.2 lc rgb 'blue', "
            "'-' with points pt 7 ps 0.5 lc rgb 'red', "
            "'-' with points pt 1 ps 1.0 lc rgb 'black', "
            "'-' with lines linestyle 1, "
            "'-' with points pt 7 ps 1.0 lc rgb 'blue'\n";
    
        m_gp << plot_str;

        // Send data
        m_gp.send1d(boost::make_tuple(cloud_vec.x, cloud_vec.y));
        m_gp.send1d(boost::make_tuple(field_nodes_vec.x, field_nodes_vec.y));
        m_gp.send1d(boost::make_tuple(quer_x, quer_y));
        m_gp.send1d(boost::make_tuple(rect_pc.x, rect_pc.y));
        m_gp.send1d(boost::make_tuple(sup_dom_x, sup_dom_y));

        // Delay one second
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
}

// Pointcloud to GpointCloud
std::vector<K::Point_2> SupportDomain::conv_pc_to_gpc(const geom::PointCloud<double>& pc)
{
    // Initialize G pointcloud
    std::vector<K::Point_2> g_pc;

    for (size_t i = 0; i < pc.pts.size(); i++)
    {
        auto pt_i = pc.pts.at(i);

        g_pc.push_back(K::Point_2(pt_i.x, pt_i.y));
    }

    return g_pc;
}

// Rectangle structure
geom::PointCloudVec<double> SupportDomain::rectangle(double x_center, double y_center,
    double width, double height)
{
    // Pointcloud initialisation
    geom::PointCloudVec<double> pc_vec;
    
    double w2 = width / 2.0; double h2 = height / 2.0;

    // Point 1 (upper left)
    pc_vec.x.push_back(x_center - w2); pc_vec.y.push_back(y_center + h2);

    // Point 2 (lower left)
    pc_vec.x.push_back(x_center - w2); pc_vec.y.push_back(y_center - h2);

    // Point 3 (lower right)
    pc_vec.x.push_back(x_center + w2); pc_vec.y.push_back(y_center - h2);

    // Point 4 (upper right)
    pc_vec.x.push_back(x_center + w2); pc_vec.y.push_back(y_center + h2);

    // Point 1 (upper left)
    pc_vec.x.push_back(x_center - w2); pc_vec.y.push_back(y_center + h2);

    return pc_vec;
}


// Get equal index
size_t SupportDomain::get_equal_idx(const K::Point_2& inter_point, const
    std::vector<K::Point_2>& data_pts, double tol)
{
    // Initialize idx
    size_t idx = 0;

    for (size_t i = 0; i < data_pts.size(); i++)
    {
        // Difference in x and y
        double diff_x = std::abs(inter_point.hx() - data_pts.at(i).hx());
        double diff_y = std::abs(inter_point.hy() - data_pts.at(i).hy());

        if (diff_x <= tol && diff_y <= tol)
        {
            idx = i;
            // break;
        }
    }
    return idx;
}