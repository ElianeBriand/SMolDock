//
// Created by eliane on 13/11/18.
//

#ifndef SMOLDOCK_AMINOACID_H
#define SMOLDOCK_AMINOACID_H


#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <tuple>
#include <memory>


namespace SmolDock {


    class Atom;

    class AminoAcid {

        friend class Protein;
    public:

        enum class AAType {
            alanine,
            arginine,
            asparagine,
            aspartate,
            cysteine,
            glutamate,
            glutamine,
            glycine,
            histidine,
            isoleucine,
            leucine,
            lysine,
            methionine,
            phenylalanine,
            proline,
            serine,
            threonine,
            tryptophan,
            tyrosine,
            valine,
            selenocysteine,
            pyrrolysine,
            heteroatom, /* Used for convenience */
            unknown
        };

        /// Internal type <=> string conversion //////
        friend std::string resTypeToString(AminoAcid::AAType t);

        friend AminoAcid::AAType stringToResType(const std::string &shorthand_or_name);


        AminoAcid() = default;

        explicit AminoAcid(const std::string &AA3LettersShorthand);

    protected:
        static std::set<std::tuple<AminoAcid::AAType, std::string, std::string> > AAShorthandSet;


    private:
        AAType type;

        std::string shorthand;
        std::string fullName;

        std::vector<std::shared_ptr<Atom> > atoms;



    };

    std::string resTypeToString(AminoAcid::AAType t);

    AminoAcid::AAType stringToResType(const std::string &shorthand_or_name);

}


#endif //SMOLDOCK_AMINOACID_H
