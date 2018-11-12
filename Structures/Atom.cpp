//
// Created by eliane on 11/11/18.
//

#include "Atom.h"

namespace SmolDock {

    unsigned int Atom::nextAtomID = 0;


    Atom::Atom(AtomType t) : type(t) {
        AtomID = nextAtomID;
        nextAtomID++;

    }

    std::string atomTypeToString(const Atom::AtomType t) {

        if (t == Atom::AtomType::carbon)
            return "Carbon";
        if (t == Atom::AtomType::oxygen)
            return "Oxygen";
        if (t == Atom::AtomType::hydrogen)
            return "Hydrogen";
    }

    Atom::AtomType Atom::getType() {
        return type;
    }

    unsigned int Atom::getAtomID() {
        return AtomID;
    }

    std::string Atom::getTypeString() {
        return atomTypeToString(type);
    }

}