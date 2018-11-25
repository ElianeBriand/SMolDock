//
// Created by eliane on 25/11/18.
//

#include "DockingResultPrinter.h"

namespace SmolDock {


    void DockingResultPrinter::printToConsole() {

    }

    DockingResultPrinter::DockingResultPrinter(std::shared_ptr<SmolDock::DockingResult> res) {
        this->res = res;
    }
}