#pragma once

#include <iostream>
#include <armadillo>
#include <vector>

#include "geom.h"

class Mesh2D
{
public:
    Mesh2D() {};
    
    // Mesh nodes indices
    std::vector<std::size_t> node_indices;

    // Coordinates of mesh nodes
    geom::PointCloud<double> node_coords;

    // Volume elements
    std::vector<geom::Triangular2DElement> volume_elements;

    // Boundary elements
    std::vector<geom::LineElement> bound_elements;

    // Volume elements centroids
    geom::PointCloud<double> get_volume_elements_centroids(void);

    // Surface elements centroids
    geom::PointCloud<double> get_surface_elements_centroids(void);

private:

    // Volume elements centroids
    geom::PointCloud<double> m_volume_elements_centroids;

    // Surface elements centroids
    geom::PointCloud<double> m_surface_elements_centroids;
};


