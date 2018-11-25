A small protein-ligand docking software.

I'm reimplementing something similar to Autodock Vina because I find the source code to be somewhat
difficult to understand.

# Build

Use cmake

# Work with it

See main.cpp for lack of a proper interface yet...

# Licencing
SmolDock is licenced under GNU GPL version 3 or later.



It includes works from :

- ESBTL (released under GNU GPLv3)


# Dependencies

## RDKit (linking with dynamic library)

TODO : Add explicit instruction on how to build. For now, review the cmakelists.txt 
to get an idea. It's a bit weird TBH.

## Vc

Edit cmakelists.txt with your Vc install dir