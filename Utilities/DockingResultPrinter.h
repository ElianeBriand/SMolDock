//
// Created by eliane on 25/11/18.
//

#ifndef SMOLDOCK_DOCKINGRESULTPRINTER_H
#define SMOLDOCK_DOCKINGRESULTPRINTER_H

#include <iostream>
#include <memory>

#include "../Structures/Results/DockingResult.h"

namespace SmolDock {

    class DockingResultPrinter {

    public:
        explicit DockingResultPrinter(std::shared_ptr<SmolDock::DockingResult> res);

        void printToConsole();


    private:
        std::shared_ptr<SmolDock::DockingResult> res;


    };

}


#endif //SMOLDOCK_DOCKINGRESULTPRINTER_H
