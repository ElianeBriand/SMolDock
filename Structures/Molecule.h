//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_MOLECULE_H
#define SMOLDOCK_MOLECULE_H

#include <vector>
#include "Structure.h"
#include "Atom.h"
#include "Bond.h"

namespace SmolDock {

    class Molecule : Structure {
    public:
        
    private:
        std::vector<Atom> atoms;
        std::vector<Bond> bonds;
    };

}

#endif //SMOLDOCK_MOLECULE_H
