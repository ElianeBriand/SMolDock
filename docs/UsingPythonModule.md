You can use the SmolDock program as python module. If you do so, you still benefit
from C++ performance because it is only a light wrapper around the code.

# Setting up SmolDock

To import the module :
```python
    import PySmolDock as sd
```
You need to have PySmolDock.so in the current directory; or else something like this
 
```python
    import os
    import sys
    # Use  os.getcwd() + "/folder1/folder2/" for relative path
    sys.path.insert(0,"/path/to/folder/containingTheSOFile/")
```

# Importing receptors/protein

```python
# We can just import PDB file (Post processor is optional)
receptor = sd.Protein()
vinaPP = sd.getVinaPostProcessorsVector()  # Enhanced similarity to Vina behaviour
receptor.populateFromPDB("./DockingTests/COX2_Ibuprofen/4PH9_COX2_without_Ibuprofen.pdb", vinaPP)
```

## Input post processors

These objects allow the exact reproduction of specific behaviour of various program. For example, the Vina post processors
make it so cysteine sulfur are not considered hydrogen bond acceptor, a behaviour of the vina program (whereas SmolDock consider it so by default).

If you need to closely emulate atom typing behaviour, you may want to use/write input post processors.

# Importing ligand

```python
# For ligand, we have many option, like populateFromSMILES.
# This shows importing from PDB
# (due to limitation of the PDB format, smiles is needed for corrent bond order)
# As above, post processor (vinaPP) is not mandatory
mol1 = sd.Molecule()
mol1.populateFromPDB("./DockingTests/COX2_Ibuprofen/VINA_Cox2_BestRes.pdb",
                     "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O",  # SMILES hint for bond order
                     120,  # Seed for random number generator (important to keep for reproducibility)
                     vinaPP)
```

We can also ```mol1.populateFromSMILES(string smiles)``` but we then rely on RDKit conformer generator for starting 3D conformation.
This is generally fine. Importing from PDB allow you to specify exact starting 3D coordinate.

You can also use other file formats : 
```python
mol1.populateFromMol2("../DockingTests/COX2_Ibuprofen/VINA_Cox2_BestRes_Charged.mol2")
mol1.populateFromMol("../DockingTests/COX2_Ibuprofen/VINA_Cox2_BestRes_Charged.mol")

```



# Setting up docking, and running it

```python
dbSetting = sd.DockingBoxSetting()
dbSetting.type = sd.DockingBoxType.centeredAround
dbSetting.center = (10.0, 22.0, 25.0) # Center of the sphere
dbSetting.radius = 10.0 # radius of the sphere

cdengine = sd.Engine.ConformerRigidDockingEngine(5, # Number of conformer
                                                 5, # Number of "retry" per conformer
                                                 receptor,
                                                 mol1,
                                                 sd.ScoringFunctionType.Vina,
                                                 sd.GlobalHeuristicType.SimulatedAnnealing,
                                                 sd.LocalOptimizerType.L_BFGS,
                                                 1244)


cdengine.setDockingBox(dbSetting)
cdengine.setupDockingEngine()
cdengine.runDockingEngine()

dockRes = cdengine.getDockingResult()

pdbwriter = sd.PDBWriter()


for mol in dockRes.ligandPoses:
    pdbwriter.addLigand(mol)

pdbwriter.writePDB("res_py.pdb")
```
# Analysing results
Exporting to PDB :
```python
dockRes = cdengine.getDockingResult()

pdbwriter = sd.PDBWriter()


for mol in dockRes.ligandPoses:
    pdbwriter.addLigand(mol)

pdbwriter.writePDB("res_py.pdb")
```

# Export ligand as RDKit molecule

```python
# Say mol1 is a docking result
rdmol1 = mol1.getRDKitMol()
```
Then use rdmol1 as any other RDKit molecule.

If mol1 is a docking result, it contains a conformer (the last one) which is the ligand final pose.


# Troubleshooting the python module

## Duplicate to_python converters
```
_frozen_importlib:219: RuntimeWarning: to-Python converter for std::vector<double, std::allocator<double> > already registered; second conversion method ignored.
```
Mostly you can ignore this. It may occur when loading other boost_python module compiled with the same boost_python version,
wich allow interoperability (SmolDock can create python RDKit object from C++ RDKit object) but leads to conflicts as bindings
for common STL containers are contained in both python module.

## No python class registered for ...
```
TypeError: No Python class registered for C++ class RDKit::ROMol
```
You are using a function to export a SmolDock Molecule (or other structure) to an object defined by another
module (mostly, RDKit), like `mol.getRDKitMol()`. This is possible if both SmolDock and the module are compiled with the exact same boost_python
version (major, minor and subrelease version number). If it is not the case, automatic conversion will fail with the aforementioned
error. Recompile either SmolDock or RDKit with the appropriate, same boost_python version. Or do not use the RDKit export features,
the rest of SmolDock (and RDKit) will work perfectly fine : they just won't communicate. You can also do LD_PRELOAD and/or rename
.so file. This mostly works.

