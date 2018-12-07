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
        auto it = std::find_if(Atom::AtomTypeLabel.begin(), Atom::AtomTypeLabel.end(),
                               [&](const std::tuple<Atom::AtomType, std::string, std::string> &e) {
                                   return std::get<0>(e) == t;
                               });
        if (it != Atom::AtomTypeLabel.end()) {
            return std::get<1>(*it);
        }
    }

    Atom::AtomType stringToAtomType(const std::string &symbol_or_name) {
        auto it = std::find_if(Atom::AtomTypeLabel.begin(), Atom::AtomTypeLabel.end(),
                               [&](const std::tuple<Atom::AtomType, std::string, std::string> &e) {
                                   return (std::get<1>(e) == symbol_or_name) || (std::get<2>(e) == symbol_or_name);
                               });
        if (it != Atom::AtomTypeLabel.end()) {
            return std::get<0>(*it);
        } else {

            std::cout << "[!] Encountered unknown atom symbol_or_name : " << symbol_or_name << std::endl;
            return Atom::AtomType::unknown;
        }
    }


    Atom::AtomType Atom::getAtomType() {
        return type;
    }

    unsigned int Atom::getAtomID() {
        return AtomID;
    }

    std::string Atom::getTypeAsString() {
        return atomTypeToString(type);
    }

    /*
     *             hydrogen,
            carbon,
            oxygen,
            nitrogen*/

    std::set<std::tuple<Atom::AtomType, std::string, std::string> > Atom::AtomTypeLabel = {
            {Atom::AtomType::unknown,  "unknown",  "?"},
            {Atom::AtomType::hydrogen, "hydrogen", "H"},
            {Atom::AtomType::carbon,   "carbon",   "C"},
            {Atom::AtomType::oxygen,   "oxygen",   "O"},
            {Atom::AtomType::nitrogen, "nitrogen", "N"},
            {Atom::AtomType::sulfur,   "sulfur",   "S"},
            {Atom::AtomType::chlorine, "chlorine", "CL"}
    };


    Atom::Atom(const std::string &symbol_or_name, bool PDBFormat, AminoAcid::AAType resType) {
        if (PDBFormat) {
            auto first_letter = symbol_or_name.substr(0, 1); // One letter code for atom is always same as usual
            this->type = stringToAtomType(first_letter);
            this->atomClassInResidue = symbol_or_name;
            fromResidue = true;
        } else {
            this->type = stringToAtomType(symbol_or_name);
        }
        this->residueType = resType; // Default : heteroatom
        this->AtomID = nextAtomID;
        nextAtomID++;
    }

    Atom::Atom(const std::string &symbol_or_name, unsigned int id) {
        this->type = stringToAtomType(symbol_or_name);
        this->AtomID = id;
    }

    std::weak_ptr<AminoAcid> Atom::getOwningAA() {
        return this->owningAA;
    }

    void Atom::setOwningAA(std::shared_ptr<AminoAcid> &aa) {
        this->owningAA = aa;
    }

    std::tuple<double, double, double> Atom::getAtomPosition() {
        return std::make_tuple(this->x, this->y, this->z);
    }

    void Atom::setAtomPosition(std::tuple<double, double, double> pos) {
        this->x = std::get<0>(pos);
        this->y = std::get<1>(pos);
        this->z = std::get<2>(pos);
    }

    iAtom Atom::generateiAtom() {
        iAtom ret{}; // R.V.O.
        emplaceiAtom(ret);
        return ret;
    }

    void Atom::emplaceiAtom(iAtom &atom) {
        atom.atomicNum = static_cast<unsigned char>(this->type); // The enum is defined such that the value are atomic numbers
        atom.x = this->x;
        atom.y = this->y;
        atom.z = this->z;
    }

    void Atom::setAtomType(Atom::AtomType t) {
        this->type = t;
    }

    Atom::AtomVariant Atom::getAtomVariant()
    {
        return this->variant;
    }

    void Atom::setAtomVariant(AtomVariant v)
    {
        this->variant = v;
    }


    void Atom::setAtomID(unsigned int id)
    {
        this->AtomID = id;
    }



}