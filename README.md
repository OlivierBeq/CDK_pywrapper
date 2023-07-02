[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# CDK Python wrapper

Python wrapper to ease the calculation of [CDK](https://cdk.github.io/) molecular descriptors and fingerprints.

## Installation

From source:

    git clone https://github.com/OlivierBeq/CDK_pywrapper.git
    pip install ./CDK_pywrapper

with pip:

```bash
pip install CDK-pywrapper
```

### Get started

```python
from CDK_pywrapper import CDK
from rdkit import Chem
from rdkit.Chem import AllChem

smiles_list = [
  # erlotinib
  "n1cnc(c2cc(c(cc12)OCCOC)OCCOC)Nc1cc(ccc1)C#C",
  # midecamycin
  "CCC(=O)O[C@@H]1CC(=O)O[C@@H](C/C=C/C=C/[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]1OC)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)OC(=O)CC)(C)O)N(C)C)O)CC=O)C)O)C",
  # selenofolate
  "C1=CC(=CC=C1C(=O)NC(CCC(=O)OCC[Se]C#N)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N",
  # cisplatin
  "N.N.Cl[Pt]Cl"
]
mols = [Chem.AddHs(Chem.MolFromSmiles(smiles)) for smiles in smiles_list]
for mol in mols:
    _ = AllChem.EmbedMolecule(mol)

cdk = CDK()
print(cdk.calculate(mols))
```

The above calculates 287 molecular descriptors (23 1D, 200 2D and 64 3D).<br/>
:warning: Molecules are required to have conformers for descriptors to be calculated.<br/>

To obtain molecular fingerprint, one can used the following:

```python
from CDK_pywrapper import CDK, FPType
cdk = CDK(FPType.PubchemFP)
print(cdk.calculate(mols))
```

The following fingerprints can be calculated:

| FPType    | Fingerprint name                                                                   |
|-----------|------------------------------------------------------------------------------------|
| FP        | CDK fingerprint                                                                    |
| ExtFP     | Extended CDK fingerprint (includes 25 bits for ring features and isotopic masses)  |
| GraphFP   | CDK fingerprinter ignoring bond orders                                             |
| MACCSFP   | Public MACCS fingerprint                                                           |
| PubchemFP | PubChem substructure fingerprint                                                   |
| SubFP     | Fingerprint describing 307 substructures                                           |
| KRFP      | Klekota-Roth fingerprint                                                           |
| AP2DFP    | Atom pair 2D fingerprint as implemented in PaDEL                                   |
| HybridFP  | CDK fingerprint ignoring aromaticity                                               |
| LingoFP   | LINGO fingerprint                                                                  |
| SPFP      | Fingerprint based on the shortest paths between two atoms                          |
| SigFP     | Signature fingerprint                                                              |
| CircFP    | Circular fingerprint                                                               |

## Documentation

```python
class CDK(fingerprint=None, nbits=1024, depth=6):
```

Constructor of a CDK calculator for molecular descriptors or fingerprints

Parameters:

- ***fingerprint  : FPType***  
  Type of fingerprint to calculate (default: None). If None, calculate descriptors
- ***nbits  : int***  
  Number of bits in the fingerprint.
- ***depth  : int***  
  Depth of the fingerprint.
<br/>
<br/>
```python
def calculate(mols, show_banner=True, njobs=1, chunksize=1000):
```

Default method to calculate BlueDesc fingerprints.

Parameters:

- ***mols  : Iterable[Chem.Mol]***  
  RDKit molecule objects for which to obtain CDK descriptors.
- ***show_banner  : bool***  
  Displays default notice about BlueDesc.
- ***njobs  : int***  
  Maximum number of simultaneous processes.
- ***chunksize  : int***  
  Maximum number of molecules each process is charged of.
