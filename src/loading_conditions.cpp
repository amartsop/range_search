#include "./loading_conditions.h"

// External force function
arma::dvec LoadingConditions::external_force_function(const arma::dvec& x,
    const arma::dvec& q_bar)
{
    arma::dvec fb = {0.0, 0.0};

    return fb;
}
    
// External traction function
arma::dvec LoadingConditions::external_traction_function(const arma::dvec& x, const
    arma::dvec& ej)
{
    // Initialize traction
    arma::dvec to = {0.0, 0.0};
    
    double d = 12.0; 
    double p = -5.0e5;

    double i_ = (1.0 / 12.0) * pow(d, 3.0);
    
    if (x(0) >= 47.8)
    {
        to(1) = 0.5 * (p / i_) * ( pow(d, 2.0) / 4.0 - pow(x(1), 2.0) );
    }

    return to;
}

