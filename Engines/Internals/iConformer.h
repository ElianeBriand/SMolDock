//
// Created by eliane on 29/11/18.
//

#ifndef SMOLDOCK_ICONFORMER_H
#define SMOLDOCK_ICONFORMER_H

#include <vector>
#include <memory>

#include "iAtom.h"

namespace SmolDock {

    struct iConformer {
        std::unique_ptr<std::vector<iAtom>> atoms_vect;
    };

}


#endif //SMOLDOCK_ICONFORMER_H
