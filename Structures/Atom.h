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
 * SmolDock is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SmolDock.  If not, see <https://www.gnu.org/licenses/>.
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

    enum class PDBResidueVariantAssignationType;

    /*!
     * \brief Class representing "rich" atom in Molecule and Protein.
     *
     */
    class Atom {
    public:

        friend Bond;
        friend MoleculeTraversal;

        friend class Molecule;

        friend class Protein;

        friend void assignApolarCarbonFlag(std::vector<std::shared_ptr<Atom> > &atomVect);

        friend void
        assignVariantFlagsForResidueAtom(AminoAcid &residue, PDBResidueVariantAssignationType assignation_type);


        // Make the value of the enum be the atomic number
        enum class AtomType : unsigned char {
            unknown = 0,
            hydrogen = 1,

            boron = 5,
            carbon = 6,
            nitrogen = 7,
            oxygen = 8,
            fluorine = 9,

            magnesium = 12,

            silicon = 14,
            phosporus = 15,
            sulfur = 16,
            chlorine = 17,

            calcium = 20,

            manganese = 25,
            iron = 26,
            cobalt = 27,

            bromine = 35,


            iodine = 53
        };

        static_assert(sizeof(unsigned int) >= 4, "unsigned int type need to be at least 32 bits for atom flag purpose");
        enum class AtomVariant : unsigned int {
            none = 0,
            apolar = 1
                    << 0, // Hydrophobic flag (notably, identifies hydrophobic carbon vs partial-charge-carrying carbon for scoring)
            hydrogenDonor = 1 << 1,
            hydrogenAcceptor = 1 << 2,

            covalentReversibleAcceptor = 1 << 3,

            // next flag : xxx = 1 << 4;
            customAtomFlag1 = 1 << 29,
            customAtomFlag2 = 1 << 29,
            customAtomFlag3 = 1 << 30

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


        static std::set<std::tuple<Atom::AtomType, std::string, std::string, double> > AtomTypeLabel;

    protected:
        // Bonds involving this atom
        std::vector<std::weak_ptr<Bond> > bonds;

    private:
        AtomType type;
        AtomVariant variant = Atom::AtomVariant::none;
        int charge = 0;

        bool fromResidue = false;
        AminoAcid::AAType residueType = AminoAcid::AAType::heteroatom;
        std::string atomClassInResidue; // The full atom name given in the PDB (CA, CB, CZ ...)
        std::weak_ptr<AminoAcid> owningAA;

        static unsigned int nextAtomID;
        unsigned int AtomID;

        double x, y, z;
        double atomicRadius;


    };

    std::string atomTypeToString(Atom::AtomType t);

    std::string atomTypeToSymbolString(Atom::AtomType t);

    Atom::AtomType stringToAtomType(const std::string &symbol_or_name);

    double atomTypeToAtomicRadius(Atom::AtomType t);


    inline constexpr Atom::AtomVariant operator|(Atom::AtomVariant a, Atom::AtomVariant b) {
        return static_cast<Atom::AtomVariant>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
    }

    inline constexpr bool operator&(Atom::AtomVariant a, Atom::AtomVariant b) {
        return static_cast<bool>(static_cast<unsigned int>(a) & static_cast<unsigned int>(b));
    }

    inline std::string atomVariantToString(Atom::AtomVariant a) {
        std::string varStr;

        if(a & Atom::AtomVariant::apolar)
        {
            varStr.append("apolar ");
        }
        if(a & Atom::AtomVariant::hydrogenDonor)
        {
            varStr.append("hydrogenDonor ");
        }
        if(a & Atom::AtomVariant::hydrogenAcceptor)
        {
            varStr.append("hydrogenAcceptor ");
        }
        if(a & Atom::AtomVariant::covalentReversibleAcceptor)
        {
            varStr.append("covalentReversibleAcceptor ");
        }
        if(a & Atom::AtomVariant::customAtomFlag1)
        {
            varStr.append("customAtomFlag1 ");
        }
        if(a & Atom::AtomVariant::customAtomFlag2)
        {
            varStr.append("customAtomFlag2 ");
        }
        if(a & Atom::AtomVariant::customAtomFlag3)
        {
            varStr.append("customAtomFlag3 ");
        }
        return varStr;
    }

}

#endif //SMOLDOCK_ATOM_H
