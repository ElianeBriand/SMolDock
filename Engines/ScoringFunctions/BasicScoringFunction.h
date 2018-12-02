//
// Created by eliane on 02/12/18.
//

#ifndef SMOLDOCK_BASICSCORINGFUNCTION_H
#define SMOLDOCK_BASICSCORINGFUNCTION_H


#include <Engines/Internals/iConformer.h>
#include <Engines/Internals/iProtein.h>

namespace SmolDock {
    namespace Score {


        double basic_scoring_func(iConformer& conformer, iProtein& protein)
        {
            return 1.0;
        }

    }

}


#endif //SMOLDOCK_BASICSCORINGFUNCTION_H
