#include "../include/mesh2D.h"

// Volume elements centroids
geom::PointCloud<double> Mesh2D::get_volume_elements_centroids(void)
{
    // Initialize centroid coordinates
    geom::PointCloud<double> centroids;

    for (size_t i = 0; i < volume_elements.size(); i++)
    {
        // Get element i 
        auto element_i = volume_elements.at(i);

        // Initialize element centroid
        geom::Point<double> centroid_i;

        // Get nodes coordinates
        geom::Point<double> pt_node1 = node_coords.pts.at(element_i.node1_index);
        geom::Point<double> pt_node2 = node_coords.pts.at(element_i.node2_index);
        geom::Point<double> pt_node3 = node_coords.pts.at(element_i.node3_index);
        
        // Get centroid coordinates
        centroid_i.x = (pt_node1.x + pt_node2.x + pt_node3.x) / 3.0;
        centroid_i.y = (pt_node1.y + pt_node2.y + pt_node3.y) / 3.0;
        centroid_i.z = (pt_node1.z + pt_node2.z + pt_node3.z) / 3.0;

        centroids.pts.push_back(centroid_i);
    }

    return centroids;
}


// Surface elements centroids
geom::PointCloud<double> Mesh2D::get_surface_elements_centroids(void)
{
    // Initialize centroid coordinates
    geom::PointCloud<double> centroids;

    for (size_t i = 0; i < bound_elements.size(); i++)
    {
        // Get element i 
        auto element_i = bound_elements.at(i);

        // Initialize element centroid
        geom::Point<double> centroid_i;

        // Get nodes coordinates
        geom::Point<double> pt_node1 = node_coords.pts.at(element_i.node1_index);
        geom::Point<double> pt_node2 = node_coords.pts.at(element_i.node2_index);
        
        // Get centroid coordinates
        centroid_i.x = (pt_node1.x + pt_node2.x) / 2.0;
        centroid_i.y = (pt_node1.y + pt_node2.y) / 2.0;
        centroid_i.z = (pt_node1.z + pt_node2.z) / 2.0;

        centroids.pts.push_back(centroid_i);

    }

    return centroids;
}
