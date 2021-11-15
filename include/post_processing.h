#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>

#include <boost/tuple/tuple.hpp>
#include "gnuplot-iostream.h"

#include "geom.h"
#include "mesh2D.h"
#include "gp_utils.h"
#include "ansys_beam.h"

class PostProcessing 
{
    public:

        PostProcessing() {}

        // Plot mesh
        void plot_mesh(Mesh2D& mesh, const std::string& color,
            bool plot_centroids=false, bool plot_numbers=false);

    private:

        // Plot mesh geometry
        void plot_mesh_geometry(Gnuplot& gp, Mesh2D& mesh, 
            const std::string& color, bool plot_centroids=false,
            bool plot_numbers=false);
        
};