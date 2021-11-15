#pragma once 

#include <iostream>
#include <vector>
#include <armadillo>
#include "./geom.h"

class LoadingConditions
{
    public:

        // External force function
        static arma::dvec external_force_function(const arma::dvec& x, const
            arma::dvec& q_bar);
    
        // External traction function
        static arma::dvec external_traction_function(const arma::dvec& x, const
            arma::dvec& q_bar);
};