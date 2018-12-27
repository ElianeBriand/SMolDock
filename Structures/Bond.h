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

#ifndef SMOLDOCK_BOND_H
#define SMOLDOCK_BOND_H

#include <memory>
#include "Atom.h"

namespace SmolDock {

    /*!
     * \brief Class representing a bond in Molecule object.
     *
     */
    class Bond : public std::enable_shared_from_this<Bond> {

    public:

        //! Supported bond type
        enum BondType {
            singlebond,
            doublebond,
            triplebond,
            defaultbond,
            aromatic
        };

        //! Create a bond from two shared_ptr to Atom
        /*!
         * \param atom_a First atom ptr
         * \param atom_b Second atom ptr
        */
        Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b);

        //! Create a bond from two shared_ptr to Atom, specifying a particular BondID
        /*!
         *
         * Use of this constructor loses the garantee that BondID are unique across a Molecule. (and, given
         * it is implemented as a static member variable, across *all* Molecule)
         *
         * \param atom_a First atom ptr
         * \param atom_b Second atom ptr
         * \param id BondID
        */
        Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b, unsigned int id);

        //! Get the integer ID of the bond
        /*!
         * \return BondID
        */
        unsigned int getBondID();

        //! Get first end of bond
        /*!
         * \return std::shared_ptr to the first atom in the bond
         * \sa getEndB()
        */
        std::shared_ptr<Atom> getEndA();

        //! Get second end of bond
        /*!
         * \return std::shared_ptr to the second atom in the bond
         * \sa getEndA()
        */
        std::shared_ptr<Atom> getEndB();

        //! Modify the weak_ptr of the atom involved in the bond to point to this bond
        /*!
         * This publicizing of the bond allow atom-based traversal. Uses std::enable_shared_from_this, thus
         * this member function has meaningful effects only if *this* bond is stored by a std::shared_ptr.
         */
        void publicizeToAtom();

        //! Set the type of the bond
        /*!
         *
         * Exotic bond type are not supported yet (if at all), as the corresponding logic for them to
         * have meaningfull effect is not present.
         *
         * \param t Bond type
         * \sa getBondType()
        */
        void setBondType(Bond::BondType t);

        //! Get the bond type
        /*!
         * \return The bond type
         * \sa setBondType()
        */
        BondType getBondType() const;

    private:
        std::shared_ptr<Atom> bond_end_a;
        std::shared_ptr<Atom> bond_end_b;

        static unsigned int nextBondID;
        unsigned int BondID;

        BondType bondtype;

    };

}

#endif //SMOLDOCK_BOND_H
