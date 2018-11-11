//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_ATOM_H
#define SMOLDOCK_ATOM_H

namespace SmolDock {

    enum AtomType {
        carbon
    };

    class Atom {
    public:
        explicit Atom(AtomType t);

    private:
        AtomType type;
    };

}

#endif //SMOLDOCK_ATOM_H
