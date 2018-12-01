//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_ATOM_H
#define SMOLDOCK_ATOM_H

#include <vector>
#include <memory>
#include <set>
#include <algorithm>
#include <iostream>

#include "AminoAcid.h"
#include "../Engines/Internals/iAtom.h"


namespace SmolDock {

    class Bond;

    class MoleculeTraversal;

    class Atom;


    class MoleculeTraversal;


    class Atom {
    public:

        friend Bond;
        friend MoleculeTraversal;

        friend class Protein;

        // Make the value of the enum be the atomic number
        enum class AtomType : unsigned char {
            unknown = 0,
            hydrogen = 1,
            carbon = 6,
            oxygen = 8,
            nitrogen = 7,
            sulfur = 16,
            chlorine = 17,
        };

        /** Too many cases, for now we just store as-is in a string
        enum class AtomClassInResidue {
            unknown,
            C, O, N,
            CA, CB, CG,
            CZ, CZ2, CZ3,
            CH2,
            CD1, CD2,
            CE1, CE2, CE3,
            CG2,
            NE, NE1, NE2,
            NH1, NH2,
            NZ,
            OE1,
            OD1, OD2,
            OG1,
            SG
        };
        */


        /// Internal type <=> string conversion //////
        friend std::string atomTypeToString(Atom::AtomType t);

        friend Atom::AtomType stringToAtomType(const std::string &symbol_or_name);

        /// Constructors //////

        explicit Atom(AtomType t);

        explicit Atom(const std::string &symbol_or_name, bool PDBFormat = false,
                      AminoAcid::AAType resType = AminoAcid::AAType::heteroatom);

        /* If you use these constructors, unique AtomID is not guaranteed */
        Atom(AtomType t, unsigned int id);

        Atom(const std::string &symbol_or_name, unsigned int id);


        AtomType getType();

        std::string getTypeString();

        unsigned int getAtomID();

        std::weak_ptr<AminoAcid> getOwningAA();

        void setOwningAA(std::shared_ptr<AminoAcid> &aa);

        std::tuple<double, double, double> getAtomPosition();

        void setAtomPosition(std::tuple<double, double, double> pos);

        iAtom generateiAtom();

        void emplaceiAtom(iAtom &atom);

    protected:
        // Bonds involving this atom
        std::vector<std::weak_ptr<Bond> > bonds;

    private:
        AtomType type;

        bool fromResidue = false;
        AminoAcid::AAType residueType = AminoAcid::AAType::heteroatom;
        std::string atomClassInResidue; // The full atom name given in the PDB (CA, CB, CZ ...)
        std::weak_ptr<AminoAcid> owningAA;

        static unsigned int nextAtomID;
        unsigned int AtomID;

        double x, y, z;

        static std::set<std::tuple<Atom::AtomType, std::string, std::string> > AtomTypeLabel;
    };

    std::string atomTypeToString(Atom::AtomType t);

    Atom::AtomType stringToAtomType(const std::string &symbol_or_name);

}

#endif //SMOLDOCK_ATOM_H
