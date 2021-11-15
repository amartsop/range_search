#include "../include/geometry_model.h"

GeometryModel::GeometryModel(const
    std::vector<SupportDomain::SupportDomainPoint>& sd_pts,
    const geom::RPIMParameters& rpim_params)
{
    // Initialize support domain points
    m_sd_pts = sd_pts;

    // Set support domain search parameters
    m_search_params = rpim_params;
    
    // Get elasticity matrix
    m_c_mat = Material::get_elasticity_matrix();
}

// Update strain and deformation
void GeometryModel::update(size_t idx, const arma::dvec& q_bar, double tol)
{
    // Generate the support domain "s"
    auto sup_dom_s = m_sd_pts.at(idx);

    // Get Ls mapping matrix
    m_ls_mat = global_to_local_mapping_matrix(sup_dom_s, q_bar.n_rows);

    // /* Nodal coordinates; for the support domain s */
    // // Mapping from global coordinates to local cordinates
    arma::dvec es = m_ls_mat * q_bar;

    /* Shape function; Calculate at the interest point for the support domain s */
    // Calculate support function on local support domain
    ShapeFunction sf_s(m_search_params.as, m_search_params.dc, m_search_params.q,
        m_sf_ms);

    // Get the interest point
    m_x_inter = sup_dom_s.point_coords;

    // Calculate shape function quantities
    auto shape_function_s = sf_s.calculate(m_x_inter, sup_dom_s.support_coords);

    /* Deformation; Calculate deformation at the interest point for the 
    support domain s*/
    m_deformation = shape_function_s.phis_mat * es;

    // Calculate deformation jacobian    
    m_deformation_jac = shape_function_s.phis_mat * m_ls_mat;

    /* Strain; Calculate strain and strain jacobian at the interest point for the 
    support domain s*/
    Strain strain;
    strain.set_shape_function(shape_function_s);
    strain.update(es);
    m_strain = strain.get_strain_vector();

    // Update ds matrix
    m_ds_matrix = strain.get_ds_matrix();

    // Calculate strain jacobian
    m_strain_jac = m_ds_matrix * m_ls_mat;
}

// Calculate f_el function
arma::dvec GeometryModel::f_el_function(const arma::dvec& x,
    const arma::dvec& q_bar)
{
   return - m_strain_jac.t() * m_c_mat * m_strain;
}

// Calculate fbex function
arma::dvec GeometryModel::f_bex_function(const arma::dvec& x,
    const arma::dvec& q_bar)
{
    return m_deformation_jac.t() * LoadingConditions::external_force_function(x, q_bar);
}

// Calculate ftex function
arma::dvec GeometryModel::f_tex_function(const arma::dvec& x,
    const arma::dvec& q_bar)
{
    return m_deformation_jac.t() * LoadingConditions::external_traction_function(x, q_bar);
}

// Calculate stifness matrix
arma::dmat GeometryModel::k_el_function(const arma::dvec& x, const arma::dvec& q_bar)
{
    return m_ls_mat.t() * m_ds_matrix.t() * m_c_mat * m_ds_matrix * m_ls_mat;
}

// Global to local coordinates mapping
arma::umat GeometryModel::global_to_local_mapping_matrix(const
    SupportDomain::SupportDomainPoint& sup_dom_s, size_t q_bar_size)
{
    // Get the number of support domain points
    size_t ns = sup_dom_s.support_indices.size();

    // Initialize ls matrix
    arma::umat ls_mat = arma::zeros<arma::umat>(2 * ns, q_bar_size);

    // Initialize i2 
    arma::umat i2 = arma::eye<arma::umat>(2, 2);

    for (size_t i = 0; i < ns; i++)
    {
        // Get support node index
        auto idx = sup_dom_s.support_indices.at(i);

        // Set ls columns
        ls_mat(arma::span(2*i, 2*i+1), arma::span(2*idx, 2*idx+1)) = i2;
    }
    
    return ls_mat;
}
