//
// Created by eliane on 26/12/18.
//

#ifndef SMOLDOCK_RESCORER_H
#define SMOLDOCK_RESCORER_H

#include <Structures/Protein.h>
#include <Structures/Molecule.h>
#include <Engines/Internals/iTransform.h>
#include <Engines/ScoringFunctions/ScoringFunctionFactory.h>
namespace SmolDock {

    //! A simple wrapper to get docking scores from various score function for a given ligand-protein configuration
    class ReScorer {
    public:
        ReScorer(Protein &prot, Molecule &mol, Score::ScoringFunctionType scorFuncType);

        bool prepare();

        double getScore();

    private:
        Protein &protein;
        Molecule &molecule;
        Score::ScoringFunctionType scoringFunctionType;

        std::shared_ptr<Score::ScoringFunction> scoringFunction;

        iProtein iprotein;
        iConformer iconformer;
        iTransform itransform;

        bool prepared = false;
    };

}

#endif //SMOLDOCK_RESCORER_H
