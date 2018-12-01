#include <iostream>
// #include <gsl/gsl_sf_bessel.h>


#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
/*
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/PreprocessRxn.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
 */
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/ExplicitBitVect.h>

#include <ESBTL/default.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <ESBTL/compressed_ifstream.h>

#define SMOLDOCK_VERBOSE_DEBUG

#include "Structures/Molecule.h"
#include "Utilities/MoleculeTraversal.h"
#include "Structures/Protein.h"
#include "Utilities/DockingResultPrinter.h"
#include "Engines/ConformerRigidDockingEngine.h"

#include <Vc/Vc>


int main() {

    SmolDock::Protein prot;
    prot.populateFromPDB("1dpx.pdb");

    SmolDock::Molecule mol;
    mol.populateFromSMILES("[CH3]O[CH2][OH]");


    SmolDock::Engine::ConformerRigidDockingEngine docker(100); // Use 100 conformers for docking

    docker.setProtein(&prot);
    docker.setLigand(&mol);
    docker.setDockingBox(SmolDock::Engine::AbstractDockingEngine::DockingBoxSetting::everything);
    docker.setRandomSeed(3984);

    if (!docker.setupDockingEngine()) {
        std::cout << "Error while doing engine setup" << std::endl;
    }

    docker.runDockingEngine();

    std::shared_ptr<SmolDock::DockingResult> res = docker.getDockingResult();

    SmolDock::DockingResultPrinter printer(res);

    printer.printToConsole();

    return 0;

    /*
    SmolDock::MoleculeTraversal tr(m);

    tr.printTraversal();

    SmolDock::Molecule m2("COCO");
    SmolDock::MoleculeTraversal tr2(m2);
    tr2.printTraversal();


    auto fingerprint_mol1 = RDKit::LayeredFingerprintMol(*m.getInternalRWMol());
    auto fingerprint_mol2 = RDKit::LayeredFingerprintMol(*m2.getInternalRWMol());

    std::cout << "Test ! " << (*fingerprint_mol1 == *fingerprint_mol2) << std::endl;



    return 0;
    */
}