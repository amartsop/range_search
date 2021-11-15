#pragma once

#include <iostream>
#include <armadillo>
#include <vector>


class Material
{
    public:

        Material(){};

        static arma::dmat get_elasticity_matrix(void);

    private:
        
};