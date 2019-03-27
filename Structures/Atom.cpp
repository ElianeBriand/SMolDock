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

#include "Atom.h"

#include <boost/log/trivial.hpp>


namespace SmolDock {

    unsigned int Atom::nextAtomID = 0;


    Atom::Atom(AtomType t) : type(t) {
        this->AtomID = nextAtomID;
        this->nextAtomID++;

        this->atomicRadius = atomTypeToAtomicRadius(t);

    }


    Atom::Atom(Atom::AtomType t, unsigned int id) : type(t) {
        this->AtomID = id;

        this->atomicRadius = atomTypeToAtomicRadius(t);

    }

    std::string atomTypeToString(const Atom::AtomType t) {
        auto it = std::find_if(Atom::AtomTypeLabel.begin(), Atom::AtomTypeLabel.end(),
                               [&](const std::tuple<Atom::AtomType, std::string, std::string, double> &e) {
                                   return std::get<0>(e) == t;
                               });
        if (it != Atom::AtomTypeLabel.end()) {
            return std::get<1>(*it);
        } else {
            BOOST_LOG_TRIVIAL(error) << "Encountered unknown atom AtomType while converting to string: "
                                     << static_cast<unsigned char>(t);
            return std::string("");
        }
    }

    std::string atomTypeToSymbolString(const Atom::AtomType t) {
        auto it = std::find_if(Atom::AtomTypeLabel.begin(), Atom::AtomTypeLabel.end(),
                               [&](const std::tuple<Atom::AtomType, std::string, std::string, double> &e) {
                                   return std::get<0>(e) == t;
                               });
        if (it != Atom::AtomTypeLabel.end()) {
            return std::get<2>(*it);
        } else {
            BOOST_LOG_TRIVIAL(error) << "Encountered unknown atom AtomType while converting to symbol string: "
                                     << static_cast<unsigned char>(t);
            return std::string("");
        }
    }

    Atom::AtomType stringToAtomType(const std::string &symbol_or_name) {
        auto it = std::find_if(Atom::AtomTypeLabel.begin(), Atom::AtomTypeLabel.end(),
                               [&](const std::tuple<Atom::AtomType, std::string, std::string, double> &e) {
                                   return (std::get<1>(e) == symbol_or_name) || (std::get<2>(e) == symbol_or_name);
                               });
        if (it != Atom::AtomTypeLabel.end()) {
            return std::get<0>(*it);
        } else {
            BOOST_LOG_TRIVIAL(error) << "Encountered unknown atom symbol_or_name while converting from string: "
                                     << symbol_or_name;
            return Atom::AtomType::unknown;
        }
    }

    double atomTypeToAtomicRadius(Atom::AtomType t) {
        auto it = std::find_if(Atom::AtomTypeLabel.begin(), Atom::AtomTypeLabel.end(),
                               [&](const std::tuple<Atom::AtomType, std::string, std::string, double> &e) {
                                   return std::get<0>(e) == t;
                               });
        if (it != Atom::AtomTypeLabel.end()) {
            return std::get<3>(*it);
        } else {
            BOOST_LOG_TRIVIAL(error) << "Encountered unknown AtomType while finding atomic radius: "
                                     << static_cast<unsigned char>(t);
            return 0;
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


    // Guestimated atomic raddi from https://www.researchgate.net/figure/Atomic-radii-used-in-Vina-and-Vinardo-scoring-functions-CA-are-aromatic-carbons-Values_fig12_303027182
    // TODO : get real values for atomic radii
    std::set<std::tuple<Atom::AtomType, std::string, std::string, double> > Atom::AtomTypeLabel = {
            {Atom::AtomType::unknown,   "unknown",    "?",  0.0},
            {Atom::AtomType::hydrogen,  "hydrogen",   "H",  1.2},
            {Atom::AtomType::carbon,    "carbon",     "C",  1.9},
            {Atom::AtomType::oxygen,    "oxygen",     "O",  1.7},
            {Atom::AtomType::nitrogen,  "nitrogen",   "N",  1.8},
            {Atom::AtomType::sulfur,    "sulfur",     "S",  2.0},
            {Atom::AtomType::chlorine,  "chlorine",   "CL", 1.8},
            {Atom::AtomType::phosporus, "phosphorus", "P",  2.1},
            {Atom::AtomType::fluorine,  "fluorine",   "F",  1.5},
            {Atom::AtomType::bromine,   "bromine",    "BR",  2.0},
            {Atom::AtomType::iodine,    "iodine",     "I",  2.2},
            {Atom::AtomType::iron,      "iron",       "FE", 6.5}, // Change this
            {Atom::AtomType::cobalt,    "cobalt",     "CO", 6.5}, // Change this
            {Atom::AtomType::cobalt,    "manganese",  "MN", 7.0}, // Change this
            {Atom::AtomType::calcium,   "calcium",    "CA", 7.9}, // Change this
            {Atom::AtomType::magnesium, "magnesium",  "MG", 7.9}, // Change this
            {Atom::AtomType::silicon,   "silicon",    "SI", 1.1}, // Change this
            {Atom::AtomType::boron,     "boron",      "B",  1.6}, // Change this
    };


    Atom::Atom(const std::string &symbol_or_name, bool PDBFormat, AminoAcid::AAType resType) {
        if (PDBFormat) {
            auto first_letter = symbol_or_name.substr(0, 1); // One letter code for atom is always same as usual
            this->type = stringToAtomType(first_letter);
            this->rawPDBAtomName = symbol_or_name;
            fromResidue = true;
        } else {
            this->type = stringToAtomType(symbol_or_name);
        }
        this->residueType = resType; // Default : heteroatom
        this->AtomID = nextAtomID;
        this->atomicRadius = atomTypeToAtomicRadius(this->type);
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


    void Atom::setAtomType(Atom::AtomType t) {
        this->type = t;
    }

    Atom::AtomVariant Atom::getAtomVariant() {
        return this->variant;
    }

    void Atom::setAtomVariant(AtomVariant v) {
        this->variant = v;
    }


    void Atom::setAtomID(unsigned int id) {
        this->AtomID = id;
    }

    unsigned int Atom::getAtomicNumber() {
        return static_cast<unsigned int>(this->type);
    }

    unsigned int Atom::getAtomVariantAsUnderlyingType() {
        return static_cast<unsigned int>(this->variant);
    }

    double Atom::getAtomicRadius() {
        return this->atomicRadius;
    }

    void Atom::setAtomicRadius(double r) {
        this->atomicRadius = r;

    }

    double Atom::getCharge() const {
        return this->charge;
    }

    void Atom::setCharge(double ch) {
        this->charge = ch;
    }

    std::string Atom::getAtomSymbol() const {
        return atomTypeToSymbolString(this->type);
    }

    const std::string& Atom::getRawPDBAtomName() const{
        return this->rawPDBAtomName;
    }


}