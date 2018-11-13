#include <iostream>
// #include <gsl/gsl_sf_bessel.h>

#include "Structures/Molecule.h"
#include "Utilities/MoleculeTraversal.h"

int main() {

    SmolDock::Molecule m;
    m._dev_populateSampleMolecule();
    SmolDock::MoleculeTraversal tr(m);

    tr.printTraversal();

    return 0;
}