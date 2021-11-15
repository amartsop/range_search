#include "./post_processing.h"

// Plot mesh
void PostProcessing::plot_mesh(Mesh2D& mesh,
    const std::string& color, bool plot_centroids, bool plot_numbers)
{
    // Initialize gnuplot
    Gnuplot gp;

    // Plot mesh geometry
    plot_mesh_geometry(gp, mesh, color, plot_centroids, plot_numbers);
}


// Plot mesh geometry
void PostProcessing::plot_mesh_geometry(Gnuplot& gp, Mesh2D& mesh,
    const std::string& color, bool plot_centroids, bool plot_numbers)
{
    // Point clouds for x and y
    std::vector<double> pc_vec_x, pc_vec_y;
    gp_utils::Range<double> range_x, range_y;

    // Geometry range
    for (size_t i = 0; i < mesh.node_coords.pts.size(); i++)
    {
        pc_vec_x.push_back(mesh.node_coords.pts.at(i).x);
        pc_vec_y.push_back(mesh.node_coords.pts.at(i).y);
    }

    // Set ranges 
    range_x = gp_utils::find_vector_range(pc_vec_x);
    range_y = gp_utils::find_vector_range(pc_vec_y);

    // Initialize global vector
    std::vector<double> line_vec_x_global, line_vec_y_global;

    // Set plot properties
    std::string single_quote = R"(')";
    gp << "set style line 1 linecolor rgb " + single_quote + color +
    single_quote + " linetype 1 linewidth 1 pointtype 7 pointsize 0.5\n";
    gp << "set key noautotitle\n";

    // Generate gnuplot string 
    std::string gp_str = "plot ";
    int loop_size = mesh.volume_elements.size();

    for (size_t i = 0; i < loop_size; i++)
    {
        if ((i < loop_size - 1) || plot_centroids)
        {
            gp_str += "'-' with linespoints linestyle 1, ";
        }
        else 
        {
            gp_str += "'-' with linespoints linestyle 1\n";
        }
    }

    // Initialize plos
    gp << "set xrange " + gp_utils::plot_range(range_x) + "\n";
    gp << "set yrange " + gp_utils::plot_range(range_y) + "\n";

    if (plot_numbers)
    {
        // Loop through elements and plot them
        for (size_t i = 0; i < loop_size; i++)
        {
            // Get Element i
            auto element_i = mesh.volume_elements.at(i);
            auto node1_index = element_i.node1_index;
            auto node2_index = element_i.node2_index;
            auto node3_index = element_i.node3_index;

            // Node coords
            auto node1_coords = mesh.node_coords.pts.at(node1_index);
            auto node2_coords = mesh.node_coords.pts.at(node2_index);
            auto node3_coords = mesh.node_coords.pts.at(node3_index);

            std::string single_quote = R"(')";

            std::string label_str1 = "set label " + single_quote +
                std::to_string(node1_index) + single_quote + " at " + 
                std::to_string(node1_coords.x) + "," +
                std::to_string(node1_coords.y) + "\n";

            std::string label_str2 = "set label " + single_quote +
                std::to_string(node2_index) + single_quote + " at " + 
                std::to_string(node2_coords.x) + "," +
                std::to_string(node2_coords.y) + "\n";

            std::string label_str3 = "set label " + single_quote +
                std::to_string(node3_index) + single_quote + " at " + 
                std::to_string(node3_coords.x) + "," +
                std::to_string(node3_coords.y) + "\n";

            gp << label_str1;
            gp << label_str2;
            gp << label_str3;
        }
    }

    // Initialize centroids coordinates
    geom::PointCloudVec<double> vol_centroids_vec, surf_centroids_vec;

    std::string ha;

    if (plot_centroids)
    {
        // Get volume elements centroids
        geom::PointCloud<double> vol_centroids = mesh.get_volume_elements_centroids();

        // Get surface elements centroids
        geom::PointCloud<double> surf_centroids = mesh.get_surface_elements_centroids();
        
        // Convert to vectors
        vol_centroids_vec = geom::conv_pc_to_pc_vec(vol_centroids);

        // Initialize surface element centroids
        surf_centroids_vec = geom::conv_pc_to_pc_vec(surf_centroids);

        gp_str += "'-' with points pt 2 ps 1 lc rgb 'red', "
            "'-' with points pt 2 ps 1.0 lc rgb 'blue'\n";
    }
    
    // Send gp_str to gp handle
    gp << gp_str;

    // Loop through elements and plot them
    for (size_t i = 0; i < loop_size; i++)
    {
        // Get Element i
        auto element_i = mesh.volume_elements.at(i);
        auto node1_index = element_i.node1_index;
        auto node2_index = element_i.node2_index;
        auto node3_index = element_i.node3_index;

        // Node coords
        auto node1_coords = mesh.node_coords.pts.at(node1_index);
        auto node2_coords = mesh.node_coords.pts.at(node2_index);
        auto node3_coords = mesh.node_coords.pts.at(node3_index);

        // Line
        std::vector<double> line_vec_x{node1_coords.x, node2_coords.x,
            node3_coords.x, node1_coords.x};

        std::vector<double> line_vec_y{node1_coords.y, node2_coords.y,
            node3_coords.y, node1_coords.y};

        // Plot
        gp.send1d(boost::make_tuple(line_vec_x, line_vec_y));
    }

    if (plot_centroids)
    {
        gp.send1d(boost::make_tuple(vol_centroids_vec.x, vol_centroids_vec.y));
        gp.send1d(boost::make_tuple(surf_centroids_vec.x, surf_centroids_vec.y));
    }
}