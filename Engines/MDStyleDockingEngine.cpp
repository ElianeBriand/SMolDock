//
// Created by eliane on 25/11/18.
//

#include "MDStyleDockingEngine.h"


namespace SmolDock {


    bool MDStyleDockingEngine::setProtein(const Protein *p) {
        this->protein = p;
        return true;
    }

    bool MDStyleDockingEngine::setMolecule(const Molecule *m) {
        this->molecule = m;
        return true;
    }


    void MDStyleDockingEngine::runDockingEngine() {

    }

    std::shared_ptr<DockingResult> MDStyleDockingEngine::getDockingResult() {
        return std::make_shared<DockingResult>();
    }

    bool MDStyleDockingEngine::setupDockingEngine() {
        return true;
    }

    bool MDStyleDockingEngine::setDockingBox(AbstractDockingEngine::DockingBoxSetting setting) {
        if (setting != DockingBoxSetting::everything) {
            std::cout << "[!] DockingBoxSetting (that is not DockingBoxSetting::everything) is not yet implemented."
                      << std::endl;
            std::cout << "[ ] Running as if DockingBoxSetting::everything was passed" << std::endl;
            return false;
        }
        return true;
    }
}