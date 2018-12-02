//
// Created by eliane on 02/12/18.
//

#ifndef SMOLDOCK_IPROTEIN_H
#define SMOLDOCK_IPROTEIN_H

#include <vector>
#include <memory>

#include "iAtom.h"


namespace SmolDock {

    struct iProtein {
        std::unique_ptr<std::vector<iAtom>> atoms_vect;
    };

}


#endif //SMOLDOCK_IPROTEIN_H
