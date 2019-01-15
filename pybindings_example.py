## Remove this/modify this if you are importing the module from somewhere else
import os
import sys

sys.path.insert(0, os.getcwd() + "/cmake-build-debug")
sys.path.insert(0, "/home/eliane/Projects/rdkit-withPyBinding/RDKit_build/install/")
##



import PySmolDock as sd
from rdkit import Chem
from rdkit.Chem import AllChem

# We can just import PDB file (Post processor is optional)
receptor = sd.Protein();
vinaPP = sd.getVinaPostProcessorsVector();  # Enhanced similarity to Vina behaviour
receptor.populateFromPDB("./DockingTests/COX2_Ibuprofen/4PH9_COX2_without_Ibuprofen.pdb", vinaPP)

# For ligand, we have many option, like populateFromSMILES.
# This shows importing from PDB
# (due to limitation of the PDB format, smiles is needed for corrent bond order)
# As above, post processor (vinaPP) is not mandatory
mol1 = sd.Molecule();
mol1.populateFromPDB("./DockingTests/COX2_Ibuprofen/VINA_Cox2_BestRes.pdb",
                     "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O",  # SMILES hint for bond order
                     120,  # Seed for random number generator (important to keep for reproducibility)
                     vinaPP)


m2 = mol1.getRDKitMol()
print(m2)


cdengine = sd.Engine.ConformerRigidDockingEngine(1,
                                                 receptor,
                                                 mol1,
                                                 sd.ScoringFunctionType.VinaRigid,
                                                 sd.GlobalHeuristicType.RandomRestart,
                                                 sd.LocalOptimizerType.L_BFGS,
                                                 1244)

cdengine.setupDockingEngine();
cdengine.runDockingEngine();

dockRes = cdengine.getDockingResult();

for mol in dockRes.ligandPoses:
    print(Chem.MolToSmiles(mol.getRDKitMol()))
