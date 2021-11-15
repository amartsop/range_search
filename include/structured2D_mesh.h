#pragma once 

#include <iostream>
#include <armadillo>
#include <string>

#include "gmsh.h"
#include "geom.h"
#include "mesh2D.h"

class Structured2DMesh
{
    public:
        /* pt1 - bottom left corner,  pt2 - top right corner */
        Structured2DMesh(const geom::Point<double>& pt1,
            const geom::Point<double>& pt2, int pts_x, int pts_y);

        // Get triangular 2D mesh
        Mesh2D get_2D_mesh(void) { return m_mesh; }

        // Get nodal spacing
        double get_nodal_spacing(void) { return m_dc; }

        // Get nodal spacing x
        double get_nodal_spacing_x(void) { return m_dcx; }

        // Get nodal spacing x
        double get_nodal_spacing_y(void) { return m_dcy; }

        // Get points in x 
        int get_points_in_x(void) { return m_points_x; }

        // Get points in y
        int get_points_in_y(void) { return m_points_y; }

    private:

        // Initialize mesh 
        Mesh2D m_mesh;

        // 2D Triangle id and number of nodes per triangle (GMSH constants)
        const int m_triangle2D_id = 2;
        const int m_triangle_nodes = 3;

        // Boundary line id and number of nodes per line (GMSH consants)
        const int m_line_id = 1;
        const int m_line_nodes = 2;

       // Nodal spacing 
       double m_dc, m_dcx, m_dcy;

       // Points in x and y
       double m_points_x, m_points_y;
};

