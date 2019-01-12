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

#include <iostream>


#include "AminoAcid.h"
#include "Atom.h"


namespace SmolDock {


    AminoAcid::AminoAcid(const std::string &AA3LettersShorthand) {
        //AminoAcid::AAType::
        shorthand = AA3LettersShorthand;
        this->type = stringToResType(AA3LettersShorthand);
    }

    std::set<std::tuple<AminoAcid::AAType, std::string, std::string> > AminoAcid::AAShorthandSet = {
            {AminoAcid::AAType::heteroatom,     "heteroatom",     "HETEROATOM"},
            {AminoAcid::AAType::unknown,        "unknown",        "UNKNOWN"},
            {AminoAcid::AAType::alanine,        "alanine",        "ALA"},
            {AminoAcid::AAType::arginine,       "arginine",       "ARG"},
            {AminoAcid::AAType::asparagine,     "asparagine",     "ASN"},
            {AminoAcid::AAType::aspartate,      "aspartate",      "ASP"},
            {AminoAcid::AAType::cysteine,       "cysteine",       "CYS"},
            {AminoAcid::AAType::glutamate,      "glutamate",      "GLU"},
            {AminoAcid::AAType::glutamine,      "glutamine",      "GLN"},
            {AminoAcid::AAType::glycine,        "glycine",        "GLY"},
            {AminoAcid::AAType::histidine,      "histidine",      "HIS"},
            {AminoAcid::AAType::isoleucine,     "isoleucine",     "ILE"},
            {AminoAcid::AAType::leucine,        "leucine",        "LEU"},
            {AminoAcid::AAType::lysine,         "lysine",         "LYS"},
            {AminoAcid::AAType::methionine,     "methionine",     "MET"},
            {AminoAcid::AAType::phenylalanine,  "phenylalanine",  "PHE"},
            {AminoAcid::AAType::proline,        "proline",        "PRO"},
            {AminoAcid::AAType::serine,         "serine",         "SER"},
            {AminoAcid::AAType::threonine,      "threonine",      "THR"},
            {AminoAcid::AAType::tryptophan,     "tryptophan",     "TRP"},
            {AminoAcid::AAType::tyrosine,       "tyrosine",       "TYR"},
            {AminoAcid::AAType::valine,         "valine",         "VAL"},
            {AminoAcid::AAType::selenocysteine, "selenocysteine", "SEC"},
            {AminoAcid::AAType::pyrrolysine,    "pyrrolysine",    "PYL"}
    };


    std::string resTypeToString(AminoAcid::AAType t) {
        auto it = std::find_if(AminoAcid::AAShorthandSet.begin(), AminoAcid::AAShorthandSet.end(),
                               [&](const std::tuple<AminoAcid::AAType, std::string, std::string> &e) {
                                   return std::get<0>(e) == t;
                               });
        if (it != AminoAcid::AAShorthandSet.end()) {
            return std::get<1>(*it);
        } else {
            return "";
        }
    }

    AminoAcid::AAType stringToResType(const std::string &shorthand_or_name) {
        auto it = std::find_if(AminoAcid::AAShorthandSet.begin(), AminoAcid::AAShorthandSet.end(),
                               [&](const std::tuple<AminoAcid::AAType, std::string, std::string> &e) {
                                   return (std::get<1>(e) == shorthand_or_name) ||
                                          (std::get<2>(e) == shorthand_or_name);
                               });
        if (it != AminoAcid::AAShorthandSet.end()) {
            return std::get<0>(*it);
        } else {
#ifdef SMOLDOCK_VERBOSE_DEBUG
            std::cout << "[!] Encountered unknown residue shorthand_or_name : " << shorthand_or_name << std::endl;
#endif
            return AminoAcid::AminoAcid::AAType::unknown;
        }
    }

    AminoAcid::AminoAcid(const AminoAcid::AAType &type) {
        this->type = type;
    }

    unsigned int AminoAcid::getAAId() const {
        return this->AAId;
    }

    void AminoAcid::setAAId(unsigned int id) {
        this->AAId = id;

    }


    bool AminoAcid::filliProtein(iProtein &prot, bool skipHydrogen) {
        for (auto &atom: this->atoms) {
            if (skipHydrogen && atom->getAtomType() == Atom::AtomType::hydrogen)
                continue;

            prot.type.push_back(atom->getAtomicNumber());
            prot.variant.push_back((unsigned int) atom->getAtomVariant()); // Not implemented


            auto pos = atom->getAtomPosition();
            prot.x.push_back(std::get<0>(pos));
            prot.y.push_back(std::get<1>(pos));
            prot.z.push_back(std::get<2>(pos));

            prot.atomicRadius.push_back(atom->getAtomicRadius());

        }
        return true;
    }

    AminoAcid::AAType AminoAcid::getType() const {
        return this->type;
    }

    void AminoAcid::setType(AminoAcid::AAType t) {
        this->type = t;

    }


}