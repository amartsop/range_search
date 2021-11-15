#pragma once

#include <iostream>
#include <armadillo>
#include <vector>
#include "nanoflann.hpp"
#include "geom.h"

class KDTrees
{
    public:
        KDTrees(const geom::PointCloud<double> interest_points,
    		const geom::PointCloud<double> dataset);
        
        // // Nearest neighbour search
        // std::vector<std::vector<int>> nn_search(int nn_number,
        //     const geom::PointCloud<double> interest_points,
        //     const geom::PointCloud<double> dataset);

        // Nearest neighbour search
        std::vector<int> radius_search(int query_pt_idx, double search_radius);

    private:

        // KDPointCloud struct
        struct KDPointCloud
        {
            // Pts container
	        std::vector<geom::Point<double>> pts;

	        // Must return the number of data points
	        inline size_t kdtree_get_point_count() const { return pts.size(); }

	        // Returns the dim'th component of the idx'th point in the class:
	        // Since this is inlined and the "dim" argument is typically an immediate value, the
	        //  "if/else's" are actually solved at compile time.
	        inline double kdtree_get_pt(const size_t idx, const size_t dim) const
	        {
		        if (dim == 0) return pts[idx].x;
		        else if (dim == 1) return pts[idx].y;
		        else return pts[idx].z;
	        }

	        // Optional bounding-box computation: return false to default to a standard bbox computation loop.
	        //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	        //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	        template <class BBOX>
	        bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
        };
    
        // Support domain typedef
	    typedef nanoflann::KDTreeSingleIndexAdaptor<
		    nanoflann::L1_Adaptor<double, KDPointCloud>,
		    KDPointCloud, 3> m_kd_tree;

		// KD interest points and dataset
		KDPointCloud m_kd_interest_points, m_kd_dataset;
};