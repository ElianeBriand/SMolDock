//
// Created by briand on 3/8/19.
//

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include "RotatableBondRemover.h"

#include <boost/log/trivial.hpp>


namespace SmolDock::InputModifier {


    void RotatableBondRemover::postProcessAtomFromLigand(SmolDock::Atom &atom) {
        // We have no special ligand-related work-arounds

    }

    void RotatableBondRemover::postProcessAtomFromProtein(SmolDock::Atom &atom, SmolDock::AminoAcid &residue) {
        // We have no special protein-related work-arounds
    }

    std::vector<std::tuple<int,int>> RotatableBondRemover::deselectRotatableBonds(std::shared_ptr<RDKit::RWMol> rwmol) {
        std::vector<std::tuple<int, int>> returnValue;

        if(this->invalid == true)
            return returnValue; // SMARTS pattern is ill formed


        std::vector<RDKit::MatchVectType> matchFlagged;

        RDKit::SubstructMatch(static_cast<const RDKit::ROMol &>(*rwmol), *(this->flaggingPattern),
                              matchFlagged);

        int atomIdxFirstEnd = -1;
        int atomIdxSecondEnd = -1;
        int numMatch = 0;
        for (unsigned int i = 0; i < matchFlagged.size(); ++i) {
            for(unsigned int j = 0; j< matchFlagged[i].size(); ++j) {

                int queryIdx = std::get<0>(matchFlagged[i][j]);
                int atomIdx = std::get<1>(matchFlagged[i][j]);

                if (queryIdx == idxOfFirstAtomMappedAtom) // We flag
                {
                    atomIdxFirstEnd = atomIdx;
                    numMatch++;
                }

                if (queryIdx == idxOfSecondAtomMappedAtom) // We flag
                {
                    atomIdxSecondEnd = atomIdx;
                    numMatch++;
                }
            }
            if(numMatch == 2 && atomIdxFirstEnd != -1 && atomIdxSecondEnd != -1)
            {
                BOOST_LOG_TRIVIAL(info) << "RotatableBondRemover : Matched one bond to avoid flagging as rotatable : ("
                << atomIdxFirstEnd<< "," <<atomIdxSecondEnd <<")" ;
                returnValue.emplace_back(std::make_tuple(atomIdxFirstEnd,atomIdxSecondEnd));
                returnValue.emplace_back(std::make_tuple(atomIdxSecondEnd,atomIdxFirstEnd));
            }
            if(numMatch > 2)
            {
                BOOST_LOG_TRIVIAL(error) << "RotatableBondRemover : Matched more than 2 end of a bond ? Check your SMARTS pattern for specificity";
                BOOST_LOG_TRIVIAL(error) << "RotatableBondRemover : SMARTS = " << this->SMARTS_pattern;
                BOOST_LOG_TRIVIAL(error) << "RotatableBondRemover : Last matched index first end of bond  = " << atomIdxFirstEnd << " (-1 = not matched)";
                BOOST_LOG_TRIVIAL(error) << "RotatableBondRemover : Last matched index second end of bond = " << atomIdxSecondEnd;
            }

            numMatch = 0;
        }


        return returnValue;

    }

    RotatableBondRemover::RotatableBondRemover(std::string SMARTS_pattern_) :
    SMARTS_pattern(SMARTS_pattern_),
    flaggingPattern(RDKit::SmartsToMol(this->SMARTS_pattern)){

        RDKit::ROMol::VERTEX_ITER it , end;
        boost::tie( it , end ) = this->flaggingPattern->getVertices();
        while( it != end ) {
            const RDKit::Atom* atom = (*(this->flaggingPattern))[*it];
            int atomMapNum = atom->getAtomMapNum();
            if(atomMapNum == 1)
            {
                idxOfFirstAtomMappedAtom = *it;
            }
            if(atomMapNum == 2)
            {
                idxOfSecondAtomMappedAtom = *it;
            }
            ++it;
        }

        if(idxOfFirstAtomMappedAtom == -1 || idxOfSecondAtomMappedAtom != -1)
        {
            BOOST_LOG_TRIVIAL(error) << "Attempting to use RotatableBondRemover with pattern : " << this->SMARTS_pattern;
            BOOST_LOG_TRIVIAL(error) << "Which does not contain two adjacent mapped atoms numbered 1 and 2 (eg COC[C:1][C:2])";
            BOOST_LOG_TRIVIAL(error) << "No rotatable bonds will be removed by this instance.";
            this->invalid = true;
        }

    }
}