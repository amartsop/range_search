#include "../include/model.h"


// Set the integration mesh
void Model::set_integration_mesh(const Geometry2DMesh::Triangular2DMesh& mesh)
{
    // Assign integration mesh to member variable
    m_integr_mesh = mesh;
}

// Integrate function
void Model::integrate(void)
{
    /******************* 2D Triangular Elements *****************/

    // Generate quadrature points for 2D triangle element
    geom::PointCloud<double> quadr_pts_t =
        m_gq_triangle.generate_quadrature_points(m_integr_mesh,
        GQTriangle2D<Model>::LINEAR);

    // Evaluate integral
    double surf_int = m_gq_triangle.evaluate(*this, &Model::f);

    /******************* 1D Line Elements *****************/

    // Generate quadrature points for 1D line element
    geom::PointCloud<double> quadr_pts_l =
        m_gq_line.generate_quadrature_points(m_integr_mesh, 4);

    // Evaluate integral
    double line_int = m_gq_line.evaluate(*this, &Model::f);

    std::cout << "surface integral: " << surf_int << std::endl;
    std::cout << "line integral: " << line_int << std::endl;
}


// Integrand function
double Model::f(const arma::dvec& x)
{
    // return pow(x(0), 2.0) * x(1) + 3;
    return 1.0;
}

