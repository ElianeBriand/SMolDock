#include <iostream>
#include <gsl/gsl_sf_bessel.h>

#include "Structures/Molecule.h"
#include "Utilities/MoleculeTraversal.h"

int main() {
    double x = 5.0;
    double y = gsl_sf_bessel_J0 (x);
    printf ("J0(%g) = %.18e\n", x, y);

    SmolDock::Molecule m;
    m._dev_populateSampleMolecule();
    SmolDock::MoleculeTraversal tr(m);

    tr.printTraversal();

    return 0;
}