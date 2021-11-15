#pragma once 

#include <iostream>
#include <vector>

namespace geom {

    // Point struct
    template <typename T>
    struct Point 
    {
        T x, y, z;
    };

    template<typename T>
    struct PointCloudVec
    {
        std::vector<T> x, y, z;
    };

    // PointCloud struct
    template <typename T>
    struct PointCloud
    {
        // Pts container
	    std::vector<geom::Point<T>>  pts;
    };


    // Support domain struct
    struct SupportDomain
    {
        double xc, yc, zc, r;
    };

    // Convert PointCloud to PointCloudVec
    template <typename T>
    PointCloudVec<T> conv_pc_to_pc_vec(const PointCloud<T>& pc)
    {
        // Initialize point cloud vec
        PointCloudVec<T> pc_vec;

        for(size_t i = 0; i < pc.pts.size(); i++)
        {
           pc_vec.x.push_back(pc.pts.at(i).x) ;
           pc_vec.y.push_back(pc.pts.at(i).y) ;
           pc_vec.z.push_back(pc.pts.at(i).z) ;
        }

        return pc_vec;
    }

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

    // RPIM parameters 
    struct RPIMParameters{

        // Exponential constant
        double q;

        // Nodal spacing dc
        double dc;

        // Nodal spacing dc_x
        double dc_x;

        // Nodal spacing dc_y
        double dc_y;

        // Support domain constant
        double as;
    };

    // Get indices of sorted array
    template <typename T>
    std::vector<size_t> sorted_indices(const std::vector<T> &v)
    {
        // initialize original index locations
        std::vector<size_t> idx(v.size());
        iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        // using std::stable_sort instead of std::sort
        // to avoid unnecessary index re-orderings
        // when v contains elements of equal values 
        std::stable_sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

        return idx;
    }
}