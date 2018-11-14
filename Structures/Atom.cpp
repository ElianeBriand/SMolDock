//
// Created by eliane on 11/11/18.
//

#include "Atom.h"

namespace SmolDock {

    unsigned int Atom::nextAtomID = 0;


    Atom::Atom(AtomType t) : type(t) {
        this->AtomID = nextAtomID;
        nextAtomID++;

    }


    Atom::Atom(Atom::AtomType t, unsigned int id) : type(t) {
        this->AtomID = id;
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