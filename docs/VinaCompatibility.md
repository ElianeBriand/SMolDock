Answers to question on whether the Vina compatible mode is suitable for your project

* The Vina-like family of scoring function (vina_like_*) are more than adequate to 
replace the intermolecular component of the Vina scoring function. The resulting score
is not binary-identical, but the discrepency is on the order of 10^-9 or better.
    * The intermolecular component of the score are the gauss1, gauss2, repulsion, hydrophobic, hydrogen bonding components 

* Contrary to Vina, the scoring function may be used directly to a given ligand-protein conformation : 
grid computation are not necessarily performed, and the scoring function might not use such grid.
    * We do not find this to have significant numerical effect on scores.

* Regarding the global binding energy results, SmolDock apply the formula given in the original Vina paper : 
`score = score_intermolecular / (1 + w * Nrot)`, with `Nrot` the number of rotatable bond between heavy atoms, and
`w` the weight factor (0.058459). 
    * The procedure used in Vina appears to not be exactly the same. We are working to reverse-engineer exactly
    what is it that is does.
    * In the meantime, there is a discrepency of around 10^-4 (kcal/mol) between the reported affinity score (~binding free energy)
    of both software.
    * It is currently not known, but is under investigation, whether this affects the relative ranking of ligand poses,
     or is just random noise.

* The atom types, and mecanisms of atom type assignation are differents, but do not affect
the result as neither Vina nor SmolDock has peculiar/unexpected classification of atoms
(an example of which would be idiosyncratic rule for assigning polar/apolar carbon).
    * Vina assign bonds between atom using distances (PDBQT files do not contain connectivity
     data), by comparison to typical covalent bonding distances for that atom pair.
        * The bond data is then used to assign hydrophobicity (apolar/polar carbon) and hydrogen-donor characteristics
        * Atom type, but not bonding data, is used for hydrogen-bond acceptor status
    * SmolDock uses connectivity-containing file format, or (PDB+SMILES hints) to obtain connectivity
    of ligands. For protein, bonds and atom typing are obtain from the ATOM name found in PDB file.
        * In a PDB file, atom belonging to the amino-acid polymer have ATOM record. The atom name field
        uniquely identifies an atom, when combined with the residue type, and complies with the IUPAC rule for
        amino acid nomenclature. For example, CB in a leucine, the beta carbon of a leucine, is fully caracterised
        by those two information, and its characteristics (apolar carbon, ...) can be deduced.
    * No difference is expected between these behaviour.
    
    
* More in depth analysis of the consequences of the discrepencies is underway, using larger collection of ligand, and poses
comparison.