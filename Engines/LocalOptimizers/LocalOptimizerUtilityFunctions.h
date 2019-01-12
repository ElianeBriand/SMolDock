//
// Created by eliane on 05/01/19.
//

#ifndef SMOLDOCK_LOCALOPTIMIZERUTILITYFUNCTIONS_H
#define SMOLDOCK_LOCALOPTIMIZERUTILITYFUNCTIONS_H

namespace SmolDock::Optimizer {

    inline bool hasNonZeroComponent(const arma::mat &vector) {
        for (unsigned int i = 0; i < vector.n_rows; i++) {
            if (vector(i) != 0.0)
                return true;
        }
        return false;
    }

}


#endif //SMOLDOCK_LOCALOPTIMIZERUTILITYFUNCTIONS_H
