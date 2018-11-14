//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_ATOM_H
#define SMOLDOCK_ATOM_H

#include <vector>
#include <memory>


namespace SmolDock {

    class Bond;

    class MoleculeTraversal;

    class Atom {
        friend Bond;
        friend MoleculeTraversal;

    public:
        enum AtomType {
            hydrogen,
            carbon,
            oxygen,
            nitrogen
        };

        friend std::string atomTypeToString(Atom::AtomType t);

        explicit Atom(AtomType t);

        /* If you use this constructor, unique AtomID is not guaranteed */
        Atom(AtomType t, unsigned int id);

        AtomType getType();

        std::string getTypeString();

        unsigned int getAtomID();

    protected:
        // Bonds involving this atom
        std::vector<std::weak_ptr<Bond> > bonds;

    private:
        AtomType type;

        static unsigned int nextAtomID;
        unsigned int AtomID;

    };

    std::string atomTypeToString(Atom::AtomType t);
}

#endif //SMOLDOCK_ATOM_H
