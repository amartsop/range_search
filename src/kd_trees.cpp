#include "../include/kd_trees.h"

KDTrees::KDTrees(const geom::PointCloud<double> interest_points,
    const geom::PointCloud<double> dataset)
{
    // Convert pointclouds to kd pointcloudes
    m_kd_interest_points.pts = interest_points.pts;

    // Convert dataset to kd dataset
    m_kd_dataset.pts = dataset.pts;
}

// Nearest neighbour search
std::vector<int> KDTrees::radius_search(int query_pt_idx, double search_radius)
{
    // Initialize vector of indices
    std::vector<int> indices;

    // Generate index for cloud
	m_kd_tree index(3, m_kd_dataset, nanoflann::KDTreeSingleIndexAdaptorParams(10));

    // Build index
	index.buildIndex(); 

    // Matches vector
    std::vector<std::pair<size_t, double>> ret_matches;

    // Set search parameters
    nanoflann::SearchParams params;

    // Get query point
    const double query_pt[3] = {m_kd_interest_points.pts.at(query_pt_idx).x,
        m_kd_interest_points.pts.at(query_pt_idx).y,
        m_kd_interest_points.pts.at(query_pt_idx).z};

    // Update influence domains
    const size_t nMatches = index.radiusSearch(&query_pt[0],
        search_radius, ret_matches, params);

    for (auto match : ret_matches)
    {
        indices.push_back(match.first);
    }
    
    return indices;
}




// }

// // Nearest neighbour search
// std::vector<std::vector<int>> KDTrees::nn_search( int nn_number,
//     const geom::PointCloud<double> interest_points,
//     const geom::PointCloud<double> dataset)
// {
//     // Initialize vector of indices
//     std::vector<std::vector<int>> indices;

//     // Convert pointclouds to kd pointcloudes
//     KDPointCloud kd_interest_points;
//     kd_interest_points.pts = interest_points.pts;

//     KDPointCloud kd_dataset;
//     kd_dataset.pts = dataset.pts;

//     // Generate index for cloud
// 	m_kd_tree index(3, kd_dataset, nanoflann::KDTreeSingleIndexAdaptorParams(10));

//     // Build index
// 	index.buildIndex(); 

//     // Loop through interest points
//     for (int i = 0; i < kd_interest_points.pts.size(); i++)
//     {
//         // Get query point
//         const double query_pt[3] = {kd_interest_points.pts.at(i).x,
//             kd_interest_points.pts.at(i).y, kd_interest_points.pts.at(i).z};
        
//         std::vector<size_t> ret_index(nn_number);
//         std::vector<double> out_dist_sqr(nn_number);

//         nn_number = index.knnSearch(&query_pt[0], nn_number,
//             &ret_index[0], &out_dist_sqr[0]);
        
//         // Initialize indices for interest point i
//         std::vector<int> interest_idx_i;
        
//         for (size_t i = 0; i < nn_number; i++)
//         {
//             interest_idx_i.push_back(ret_index[i]);
//         }

//         indices.push_back(interest_idx_i);
//     }

//     return indices;
// }
