#include <iostream>
// #include <gsl/gsl_sf_bessel.h>


#include "Structures/Molecule.h"
#include "Utilities/MoleculeTraversal.h"
#include "Structures/Protein.h"
#include "Utilities/DockingResultPrinter.h"
#include "Engines/ConformerRigidDockingEngine.h"



int main() {

    SmolDock::Protein prot;
    // prot.populateFromPDB("1dpx.pdb"); // Lysozyme
    prot.populateFromPDB("../DockingTests/COX2_Ibuprofen/3LN1_NoHeme_NoLigand.pdb"); // COX-2

    SmolDock::Molecule mol;
    mol.populateFromSMILES("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"); // Ibuprofen


    SmolDock::Engine::ConformerRigidDockingEngine docker(100); // Use 100 conformers for docking

    docker.setProtein(&prot);
    docker.setLigand(&mol);
    docker.setDockingBox(SmolDock::Engine::AbstractDockingEngine::DockingBoxSetting::everything);
    docker.setRandomSeed(3984);

    if (!docker.setupDockingEngine()) {
        std::cout << "[!] Error while doing engine setup" << std::endl;
        return 2;
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