//
// Created by eliane on 08/01/19.
//

#ifndef SMOLDOCK_HEURISTICINTERFACE_H
#define SMOLDOCK_HEURISTICINTERFACE_H

namespace SmolDock::Heuristics {

    class GlobalHeuristic {

    public:
        virtual bool search() = 0;

        virtual arma::mat getResultMatrix() = 0;

        virtual ~GlobalHeuristic() = default;

    private:
    };


}


#endif //SMOLDOCK_HEURISTICINTERFACE_H
