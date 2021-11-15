#include "../include/shape_function.h"


ShapeFunction::ShapeFunction(double ac, double dc, double q, int ms)
{
    // Get shape function constants
    m_ac = ac; m_dc = dc; m_q = q; m_ms = ms;
}


// Calculate
ShapeFunction::Measures ShapeFunction::calculate(const arma::dvec& x,
    const std::vector<arma::dvec>& sup_dom)
{
    // Number of sample points
    size_t ns = sup_dom.size();

    // Initialize measures stucture
    Measures phis_str;
    phis_str.phis_vec = arma::zeros<arma::dvec>(ns, 1);
    phis_str.phis_jac = arma::zeros<arma::dmat>(ns, x.n_rows);
    phis_str.phis_mat =  arma::zeros<arma::dmat>(x.n_rows, 2 * ns);

    // ri and pi vectors and jacobians
    VecJac ri_vecjac = rbf_mq(x, sup_dom);
    VecJac pi_vecjac = polynomial2D_basis(x);

    // Calculate gs matrix
    arma::dmat gs_mat = gs_matrix(sup_dom);

    // Calculate gs inverse 
    arma::dmat gs_tilde_mat = arma::inv(gs_mat);

    // Initialize i2 mat
    arma::dmat i2 = arma::eye(x.n_rows, x.n_rows);

    // PARALLELISE
    // Loop through support domain points
    for (size_t i = 0; i < ns; i++)
    {
        // Caclulate the first sum
        double sum1 = 0;
        arma::drowvec sum1_jac = arma::zeros<arma::drowvec>(1, x.n_rows);

        for (size_t j = 0; j < ns; j++)
        {
            sum1 += ri_vecjac.vec.at(j) * gs_tilde_mat.at(j, i);
            sum1_jac += ri_vecjac.jac.row(j) * gs_tilde_mat.at(j, i);
        }
        
        // Caclulate the second sum
        double sum2 = 0.0;
        arma::drowvec sum2_jac = arma::zeros<arma::drowvec>(1, x.n_rows);
        for (size_t k = 0; k < m_ms; k++)
        {
            sum2 += pi_vecjac.vec.at(k) * gs_tilde_mat.at(ns+k, i);
            sum2_jac += pi_vecjac.jac.row(k) * gs_tilde_mat.at(ns+k, i);
        }

        // Construct vector and jacobian
        phis_str.phis_vec.at(i) = sum1 + sum2;
        phis_str.phis_jac.row(i) = sum1_jac + sum2_jac;

        // Calculate Phis matrix
        phis_str.phis_mat.cols(2*i, 2*i+1) =  phis_str.phis_vec.at(i) * i2;
    }

    return phis_str;
}

// Gs matrix
arma::dmat ShapeFunction::gs_matrix(const std::vector<arma::dvec>& sup_dom)
{
    // Number of sample points
    size_t ns = sup_dom.size();

    // Initialize Rs_tilde matrix
    arma::dmat rs_tilde  = arma::zeros(ns, ns);

    // Initialize Ps_tilde matrix
    arma::dmat ps_tilde  = arma::zeros(ns, m_ms);

    // PARALLELISE
    for (size_t i = 0; i < ns; i++)
    {
        // Get the support point i
        arma::dvec x_si = sup_dom.at(i);

        // Calculate the rbf vector of x_si
        VecJac rs_x_si_vecjac = rbf_mq(x_si, sup_dom);

        // Calculate the polynomial basis vector of x_si
        VecJac ps_x_si_vecjac = polynomial2D_basis(x_si);

        // Calculate the i row of the rs_tilde matrix
        rs_tilde.row(i) = rs_x_si_vecjac.vec.t();

        // Calculate the i row of the ps_tilde matirx
        ps_tilde.row(i) = ps_x_si_vecjac.vec.t();
    }

    // Construct gs matrix
    arma::dmat gs_rows1 = arma::join_horiz(rs_tilde, ps_tilde);
    arma::dmat gs_rows2 = arma::join_horiz(ps_tilde.t(), arma::zeros(m_ms, m_ms));

    return arma::join_vert(gs_rows1, gs_rows2);
}

// Radial basis function (multi-quadrics)
ShapeFunction::VecJac ShapeFunction::rbf_mq(const arma::dvec& x,
    const std::vector<arma::dvec>& sup_dom)
{
    // Initialzie vedjac
    VecJac vec_jac; 

    // Get support domain size
    size_t ns = sup_dom.size();

    // Initialize vector and jacobian
    vec_jac.vec = arma::zeros<arma::dvec>(ns, 1);

    vec_jac.jac = arma::zeros<arma::dmat>(ns, x.n_rows);

    // PARALLELISE
    for (size_t i = 0; i < ns; i++)
    {
        // Get support domain point i 
        arma::dvec xi = sup_dom.at(i);

        // D1 function
        double d1_i = x(0) - xi(0);

        // D2 function
        double d2_i = x(1) - xi(1);

        // D function
        double di = pow(d1_i, 2.0) + pow(d2_i, 2.0) + pow(m_ac * m_dc, 2.0);

        // Calculate r_vector component
        vec_jac.vec(i) = pow(di, m_q);

        // Calculate jacobian
        arma::dvec jac_row =  {d1_i, d2_i};
        vec_jac.jac.row(i) = 2.0 * m_q * pow(di, m_q-1) * jac_row.t();
    }
    
    return vec_jac;
}

// Polynomial function (linear basis in 2D)
ShapeFunction::VecJac ShapeFunction::polynomial2D_basis(const arma::dvec& x)
{
    // Initialzie vedjac
    VecJac vec_jac; 

    if (m_ms == 3)
    {
        vec_jac.vec = {1, x(0), x(1)};
        vec_jac.jac = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
    }
    else
    {
        arma::dvec a; vec_jac.vec = a;
        arma::dmat a_jac; vec_jac.jac = a_jac; 
    }

    return vec_jac;
}