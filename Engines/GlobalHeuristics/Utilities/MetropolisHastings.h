//
// Created by briand on 2/20/19.
//

#ifndef SMOLDOCK_METROPOLISHASTINGS_H
#define SMOLDOCK_METROPOLISHASTINGS_H


#include <random>

namespace SmolDock::Heuristics {

    bool MetropolisAccept(double oldScore, double newScore, double temperature, std::mt19937 &rndGenerator);

}

#endif //SMOLDOCK_METROPOLISHASTINGS_H
