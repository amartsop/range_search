#pragma once 

#include <iostream>
#include <vector>
#include <algorithm>

namespace gp_utils {


    // Range stucture
    template <typename T>
    struct Range
    {
        T min, max, val;
    };

    // Return min and max value of a vector
    template <typename T>
    static Range<T> find_vector_range(const std::vector<T>& val_vec)
    {
        // Calculate plot range
        auto max_val = std::max_element(std::begin(val_vec), std::end(val_vec));
        auto min_val = std::min_element(std::begin(val_vec), std::end(val_vec));

        Range<T> range{*min_val, *max_val, *max_val - *min_val};
        return range;
    }

    // Calculate plot's range
    static std::string plot_range(const std::vector<double>& val_vec)
    {
        // Calculate plot range
        auto range = find_vector_range<double>(val_vec);

        // Correct min and max vals
        double min_val_corr = (range.min - range.val / 2.0);
        double max_val_corr = (range.max + range.val / 2.0);
        
        // Create string for gnuplot
        const std::string range_val = "[" + std::to_string(min_val_corr) + ":" +
            std::to_string(max_val_corr) + "]";

        return range_val;
    }

    // Calculate plot's range
    template <typename T>
    static std::string plot_range(const Range<T>& range)
    {
        // Correct min and max vals
        double min_val_corr = (range.min - range.val / 2.0);
        double max_val_corr = (range.max + range.val / 2.0);
        
        // Create string for gnuplot
        const std::string range_val = "[" + std::to_string(min_val_corr) + ":" +
            std::to_string(max_val_corr) + "]";

        return range_val;
    }

}