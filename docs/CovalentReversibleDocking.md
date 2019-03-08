You can use the SmolDock program to run docking procedure taking into account covalent interaction
between select atoms.

# Example use cases
You have a serine protease, and covalent (or covalent-reversible) inhibitors. Non-covalent docking
gives you some results, but now you want to compare these poses to poses with covalent interaction,
to get an idea of whether the covalent poses are a natural continuation of the non-covalent poses. 
(This would suggest a two-phases model : non-covalent approach, setting the geometry for nucleophilic
attack -- and then covalent bonding).

This is in fact the motivating use case for this feature.

# How to

```C++

namespace sd = SmolDock

sd::Protein prot;
sd::Molecule mol;

// Apply the reversible acceptor tag to the ligand atom
mol.applyAtomVariant("c[C:1](=O)C", sd::Atom::AtomVariant::covalentReversibleAcceptor);

// Apply the special typing covalentReversibleSerineOH
prot.applySpecialResidueTyping(sd::AminoAcid::AAType::serine,221,sd::SpecialResidueTyping::covalentReversibleSerineOH);

// Use an appropriate scoring function : the covalentReversible variant to what you use for conventional docking
// Eg : sd::Score::ScoringFunctionType::VinaCovalentReversible
sd::Engine::ConformerRigidDockingEngine docker(5, 3,&prot,&mol,
                                                   sd::Score::ScoringFunctionType::VinaCovalentReversible,
                                                   sd::Heuristics::GlobalHeuristicType::SimulatedAnnealing,
                                                   sd::Optimizer::LocalOptimizerType::L_BFGS,
                                                   1244);

```

It is possible to set covalentReversibleDonor on the ligand, and use a special type that sets
 the acceptor type on the protein. Donor more of less means nucleophile, and acceptor means electrophile,
 but the terminology is not consistent for all interactions (which are not necessarily all 
 nucleophilic attack)
 
 # Data
 
 Type of covalent interaction available :
 
 | SpecialResidueTyping | Corresponding atom variant for ligand | Biological relevance |
 | --- | --- |--- |
 | covalentReversibleSerineOH | covalentReversibleAcceptor | Catalytic triad/serine hydrolase |
 | covalentReversibleCysteineSH | covalentReversibleAcceptor | Catalytic triad/cysteine hydrolase |
 
 Useful SMARTS match pattern :
 
 *  `C[C:1](=O)C` : for carbonyl electrophile. Need to adjust the neighbouring C to your chemical class.
 * `B` : boron (nice !).
 