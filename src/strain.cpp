#include "../include/strain.h"

// Set shape function
void Strain::set_shape_function(const ShapeFunction::Measures& sf)
{
    // Store shape function to member variable
    m_sf_s = sf;

    // Gen number of support domain points
    m_ns = sf.phis_vec.size();
}

// Update function
void Strain::update(const arma::dvec& es)
{
    // Phis jac
    arma::dmat phis_jac = m_sf_s.phis_jac;

    // Initialize Ds matrix 
    m_ds_matrix = arma::zeros<arma::dmat>(3, 2 * m_ns);

    // Calculate Ds matrix
    for (size_t i = 0; i < m_ns; i++)
    {
        double dphi_is_x1 = phis_jac.at(i, 0);
        double dphi_is_x2 = phis_jac.at(i, 1);

        arma::dmat di_s = {{dphi_is_x1, 0.0}, {0.0, dphi_is_x2}, 
            {dphi_is_x2, dphi_is_x1}};

        m_ds_matrix.cols(2*i, 2*i+1) = di_s;
    }

    m_strain_vector = m_ds_matrix * es;
}