/*
 * Copyright (c) 2018 Eliane Briand
 *
 * This file is part of SmolDock.
 *
 * SmolDock is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef SMOLDOCK_ATOM_H
#define SMOLDOCK_ATOM_H

#include <vector>
#include <memory>
#include <set>
#include <algorithm>
#include <iostream>

#include "AminoAcid.h"


namespace SmolDock {

    class Bond;

    class MoleculeTraversal;

    class Atom;


    class MoleculeTraversal;


    class Atom {
    public:

        friend Bond;
        friend MoleculeTraversal;
        friend class Molecule;
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

        enum class AtomVariant : unsigned int {
            unknown = 0,
            aromaticCarbon = 1,
            aromaticNitrogen = 2,
            sp3carbon = 3,
            sp2carbon = 4,
            sp1carbon = 5
        };

        /* Too many cases, for now we just store as-is in a string
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


        // //// Internal type <=> string conversion //////
        friend std::string atomTypeToString(Atom::AtomType t);
        friend std::string atomTypeToSymbolString(Atom::AtomType t);
        friend Atom::AtomType stringToAtomType(const std::string &symbol_or_name);

        // //// Properties matcher helper //////
        friend double atomTypeToAtomicRadius(Atom::AtomType t);

        // /// Constructors //////

        explicit Atom(AtomType t);

        //! Construct atom from symbol (C,O,N), lower-case atom name (carbon, hydrogen) or optionnally PDB file symbol (CA,NZ, ...)
        /*!
         *
         * \param symbol_or_name Either a one-letter symbol, lower-case atom name, or a PDB file symbol.
         * \param PDBFormat Parse as PDB file symbol.
         * \sa Atom(AtomType t)
        */
        explicit Atom(const std::string &symbol_or_name, bool PDBFormat = false,
                      AminoAcid::AAType resType = AminoAcid::AAType::heteroatom);

        /* If you use these constructors, unique AtomID is not guaranteed */
        Atom(AtomType t, unsigned int id);

        Atom(const std::string &symbol_or_name, unsigned int id);


        AtomType getAtomType();
        void setAtomType(AtomType t);
        std::string getTypeAsString();
        std::string getAtomSymbol();
        unsigned char getAtomicNumber();


        AtomVariant getAtomVariant();
        void setAtomVariant(AtomVariant v);
        unsigned int getAtomVariantAsUnderlyingType();


        unsigned int getAtomID();
        void setAtomID(unsigned int id);

        std::weak_ptr<AminoAcid> getOwningAA();
        void setOwningAA(std::shared_ptr<AminoAcid> &aa);

        std::tuple<double, double, double> getAtomPosition();
        void setAtomPosition(std::tuple<double, double, double> pos);

        double getAtomicRadius();
        void setAtomicRadius(double r);

        int getCharge();
        void setCharge(int ch);



    protected:
        // Bonds involving this atom
        std::vector<std::weak_ptr<Bond> > bonds;

    private:
        AtomType type;
        AtomVariant variant;
        int charge = 0;

        bool fromResidue = false;
        AminoAcid::AAType residueType = AminoAcid::AAType::heteroatom;
        std::string atomClassInResidue; // The full atom name given in the PDB (CA, CB, CZ ...)
        std::weak_ptr<AminoAcid> owningAA;

        static unsigned int nextAtomID;
        unsigned int AtomID;

        double x, y, z;
        double atomicRadius;

        static std::set<std::tuple<Atom::AtomType, std::string, std::string, double> > AtomTypeLabel;
    };

    std::string atomTypeToString(Atom::AtomType t);
    std::string atomTypeToSymbolString(Atom::AtomType t);

    Atom::AtomType stringToAtomType(const std::string &symbol_or_name);

    double atomTypeToAtomicRadius(Atom::AtomType t);

}

#endif //SMOLDOCK_ATOM_H
