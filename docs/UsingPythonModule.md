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


# Importing ligand

# Export ligand as RDKit molecule

```python
# Say mol1 is a docking result
rdmol1 = mol1.getRDKitMol()
```
Then use rdmol1 as any other RDKit molecule.

If mol1 is a docking result, it contains a conformer (the last one) which is the ligand final pose.

# Setting up docking, and running it


# Analysing results


# Troubleshooting the python module

## Duplicate to_python converters
```
_frozen_importlib:219: RuntimeWarning: to-Python converter for std::vector<double, std::allocator<double> > already registered; second conversion method ignored.
```
Mostly you can ignore this. It may occur when loading other boost_python module compiled with the same boost_python version,
wich allow interoperability (SmolDock can create python RDKit object from C++ RDKit object) but leads to conflicts as bindings
for common STL containers are contained in both python module.

# No python class registered for ...
```
TypeError: No Python class registered for C++ class RDKit::ROMol
```
You are using a function to export a SmolDock Molecule (or other structure) to an object defined by another
module (mostly, RDKit), like `mol.getRDKitMol()`. This is possible if both SmolDock and the module are compiled with the exact same boost_python
version (major, minor and subrelease version number). If it is not the case, automatic conversion will fail with the aforementioned
error. Recompile either SmolDock or RDKit with the appropriate, same boost_python version. Or do not use the RDKit export features,
the rest of SmolDock (and RDKit) will work perfectly fine : they just won't communicate.

