#pragma once

#include <iostream>
#include <armadillo>
#include <vector>
#include "./support_domain.h"
#include "./shape_function.h"
#include "./strain.h"
#include "./material.h"
#include "./loading_conditions.h"


class GeometryModel
{
public:
    GeometryModel(const std::vector<SupportDomain::SupportDomainPoint>& sd_pts,
        const geom::RPIMParameters& rpim_params);

    /**
    * Updates the geometric models. This includes the calculation of the strain, 
    * deformation and jacobians du/dq_bar and depsilon/dq_bar
    *
    * @param idx The idx of the point in which the deformation is evaluated 
    * x_idx = m_sd_pts.at(idx).point_coords. The function is evaluated to discrete 
    * x values (defined by idx). It can't be used for any x.
    * @param q_bar The global vector of deformations.
    */
    void update(size_t idx, const arma::dvec& q_bar, double tol=1.0e-5);

    // Get strain 
    arma::dvec get_strain(void) { return m_strain; }

    // Get deformation
    arma::dvec get_deformation(void) { return m_deformation; }
    
    // Get strain jacobian
    arma::dmat get_strain_jacobian(void) { return m_strain_jac; }

    // Get deformation jacobian
    arma::dmat get_deformation_jacobian(void) { return m_deformation_jac; } 

    // Get x interest
    arma::dvec get_x_interest(void) { return m_x_inter; }

    // Get Ds matrix    
    arma::dmat get_ds_matrix(void) { return m_ds_matrix; }

    // Get elasticity matrix
    arma::dmat get_elasticity_matrix(void) { return m_c_mat; }

public:
    // Calculate f_el function
    arma::dvec f_el_function(const arma::dvec& x, const arma::dvec& q_bar);

    // Calculate fbex function
    arma::dvec f_bex_function(const arma::dvec& x, const arma::dvec& q_bar);

    // Calculate ftex function
    arma::dvec f_tex_function(const arma::dvec& x, const arma::dvec& q_bar);

    // Calculate stifness matrix
    arma::dmat k_el_function(const arma::dvec& x, const arma::dvec& q_bar);

private:

    // Support domain points structure
    std::vector<SupportDomain::SupportDomainPoint> m_sd_pts;
    
    // Support domain radius 
    geom::RPIMParameters m_search_params;

private:
    /**
    * Calculates the mapping from global to local coordinates
    *
    * @param sup_dom_s Support domain structure for the "s" domain
    * @param q_bar_size The size of the global vector of deformations.
    * @return Local (to the support domain) mapping matrix
    */
    arma::umat global_to_local_mapping_matrix(const
        SupportDomain::SupportDomainPoint& sup_dom_s, size_t q_bar_size);

private:

    // // Shape function constants
    double m_sf_ms = 0.0;

    // Deformation vector
    arma::dvec m_deformation;

    // Strain vector
    arma::dvec m_strain;

    // Deformation jacobian
    arma::dmat m_deformation_jac;

    // Strain jacobian
    arma::dmat m_strain_jac;

    // x interest 
    arma::dvec m_x_inter;

    // Elasticity matrix
    arma::dmat m_c_mat;

    // Mapping matrix
    arma::umat m_ls_mat;

    // Ds matrix   
    arma::dmat m_ds_matrix;
};