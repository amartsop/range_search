#include <iostream>
#include <armadillo>
#include <vector>

#include "./include/structured2D_mesh.h"
#include "./include/geom.h"
#include "./include/mesh2D.h"
#include "./include/rpim2D.h"

int main(void)
{
    // Bottom left point 
    geom::Point<double> pt1{0.0, -6.0, 0.0};

    // Top right boint 
    geom::Point<double> pt2{48.0, 6.0, 0.0};

    // Structured mesh
    int iter = 6;
    int ptx_iter = iter; int pty_iter = iter;
    int pts_x = (2 * ptx_iter); int pts_y = (2 * pty_iter);
    Structured2DMesh str_mesh_raw(pt1, pt2, pts_x, pts_y);
    Mesh2D mesh_raw = str_mesh_raw.get_2D_mesh();

    // RPIM parameters
    geom::RPIMParameters rpim_params;
    
    rpim_params.q = 1.03; rpim_params.as = 4.0;
    rpim_params.dc_x = str_mesh_raw.get_nodal_spacing_x();
    rpim_params.dc_y = str_mesh_raw.get_nodal_spacing_y();
    rpim_params.dc = str_mesh_raw.get_nodal_spacing();

    // Model thickness
    double thickness = 1.0;

    // Generate model
    RPIM2D cantilever_beam;
    cantilever_beam.initialize(mesh_raw, thickness, rpim_params, true);

    /************************* Static linear analysis **********6***************/
    // Get total number of nodes
    size_t ndofs = cantilever_beam.get_dofs_num();

    // Generate initial conditions
    arma::dvec q0_bar = arma::zeros<arma::dvec>(ndofs, 1);
    
    // Update matrices
    cantilever_beam.update(q0_bar);   

    // std::cout << "Update complete!" << std::endl;
}