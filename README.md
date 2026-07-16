# 🧪 CDK Python Wrapper

<!-- Badges -->
<div align="center">

[![PyPI version](https://img.shields.io/pypi/v/CDK-pywrapper.svg)](https://pypi.org/project/CDK-pywrapper/)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/CDK-pywrapper.svg)](https://pypi.org/project/CDK-pywrapper/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/OlivierBeq/CDK_pywrapper/actions/workflows/tests.yml/badge.svg)](https://github.com/OlivierBeq/CDK_pywrapper/actions/workflows/tests.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
<br>
</div>

A simple and reliable Python wrapper for calculating [CDK](https://cdk.github.io/) (Chemistry Development Kit) molecular descriptors and fingerprints. This library takes care of installing a matching Java runtime, dispatching molecules to the bundled CDKdesc executable, and collecting the results into a tidy pandas DataFrame — so you can stay in RDKit/pandas-land.

## ✨ Features

- 🧬 **288 descriptors & 14 fingerprint types** — 223 1D/2D and 65 3D descriptors, plus FP, ExtFP, EStateFP, GraphFP, MACCSFP, PubchemFP, SubFP, KRFP, AP2DFP, HybridFP, LingoFP, SPFP, SigFP and CircFP fingerprints, straight from CDK.
- ☕ **Zero Java setup** — automatically downloads, caches and reuses a matching JRE on first use; nothing to install by hand.
- ⚡ **Parallel by design** — spread the work across multiple CPU cores with configurable `njobs`/`chunksize`, each worker running its own single-core-pinned JVM.
- 🧯 **Never silently misaligned** — molecules that fail, get skipped, or lack hydrogens/conformers are handled explicitly and reported, never dropped without a trace.
- 📊 **pandas-native output** — results come back as a ready-to-use `DataFrame`, one row per molecule.
- 🧾 **Optional CDK canonical SMILES** — get back the SMILES CDK itself parsed alongside the descriptor/fingerprint values.
- 🔀 **SD or SMILES interchange format** — choose whether molecules are handed to CDK as V2000 SD blocks or as plain SMILES.
- 🔍 **Rich metadata** — inspect every descriptor's name, description, type and dimensionality (1D/2D/3D) programmatically via `get_details()`.

## ✍️ Copyright and Citation Notice

Olivier J. M. Béquignon is **neither** the copyright holder of CDK **nor** responsible for it. The work carried out here concerns:

- the Python wrapper,
- the CDKdesc executable.

### Citing

If you use this wrapper in your research, please cite the original CDK publications in addition to this software package:

1. **Original CDK papers:**
   > Willighagen, E.L. et al. (2017), The Chemistry Development Kit (CDK) v2.0: atom typing, depiction, molecular formulas, and substructure searching. *Journal of Cheminformatics*, 9, 33.
   > [DOI: 10.1186/s13321-017-0220-4](https://doi.org/10.1186/s13321-017-0220-4)
   >
   > May, J.W. and Steinbeck, C. (2014), Efficient ring perception for the Chemistry Development Kit. *Journal of Cheminformatics*, 6, 3.
   > [DOI: 10.1186/1758-2946-6-3](https://doi.org/10.1186/1758-2946-6-3)
   >
   > Steinbeck, C. et al. (2006), Recent Developments of the Chemistry Development Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics. *Current Pharmaceutical Design*, 12(17), 2111-2120.
   > [DOI: 10.2174/138161206777585274](https://doi.org/10.2174/138161206777585274)
   >
   > Steinbeck, C. et al. (2003), The Chemistry Development Kit (CDK): An Open-Source Java Library for Chemo- and Bioinformatics. *Journal of Chemical Information and Computer Sciences*, 43(2), 493-500.
   > [DOI: 10.1021/ci025584y](https://doi.org/10.1021/ci025584y)

2. **This wrapper:**
   > Béquignon, O. J. M. *CDK_pywrapper: a Python wrapper for CDK molecular descriptors and fingerprints.* https://github.com/OlivierBeq/CDK_pywrapper

## 📦 Installation

```bash
pip install CDK-pywrapper
```

Or from source:

```bash
git clone https://github.com/OlivierBeq/CDK_pywrapper.git
pip install ./CDK_pywrapper
```

## 🛠️ Requirements

- Python 3.11+
- [RDKit](https://www.rdkit.org/docs/Install.html)

## 💡 Usage

### 1D and 2D descriptors

```python
from CDK_pywrapper import CDK
from rdkit import Chem

smiles_list = [
    # erlotinib
    "n1cnc(c2cc(c(cc12)OCCOC)OCCOC)Nc1cc(ccc1)C#C",
    # midecamycin
    "CCC(=O)O[C@@H]1CC(=O)O[C@@H](C/C=C/C=C/[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]1OC)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)OC(=O)CC)(C)O)N(C)C)O)CC=O)C)O)C",
    # selenofolate
    "C1=CC(=CC=C1C(=O)NC(CCC(=O)OCC[Se]C#N)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N",
    # cisplatin
    "N.N.Cl[Pt]Cl",
]
mols = [Chem.AddHs(Chem.MolFromSmiles(smiles)) for smiles in smiles_list]

cdk = CDK()
print(cdk.calculate(mols))
```

This calculates 223 molecular descriptors (23 1D and 200 2D).

Optionally, request the canonical SMILES CDK itself generated for each parsed molecule:

```python
print(cdk.calculate(mols, cdk_smiles=True))
```

### 3D descriptors

By default, the `ignore_3D` parameter is set to `True`, preventing any 3D descriptor from being calculated.

Should molecules with 3D coordinates be provided, one can turn on the additional 65 three-dimensional descriptors:

```python
from rdkit.Chem import AllChem

for mol in mols:
    _ = AllChem.EmbedMolecule(mol)

cdk3d = CDK(ignore_3D=False)
print(cdk3d.calculate(mols))
```

> ⚠️ A warning is raised if molecules lack hydrogens.
> ⚠️ An exception is raised if a 3D descriptor is requested for a conformer-less molecule.

```python
mol = Chem.MolFromSmiles('CCC')

cdk3d = CDK(ignore_3D=False)
print(cdk3d.calculate([mol]))
# ValueError: Cannot calculate the 3D descriptors of a conformer-less molecule
```

### Fingerprints

```python
from CDK_pywrapper import CDK, FPType

cdk = CDK(fingerprint=FPType.PubchemFP)
print(cdk.calculate(mols))
```

Fingerprint size and search depth can be tuned for the fingerprints that support it:

```python
cdk = CDK(fingerprint=FPType.FP, nbits=2048, depth=8)
print(cdk.calculate(mols))
```

The following fingerprints can be calculated:

| FPType    | Fingerprint name                                                                  |
|-----------|-------------------------------------------------------------------------------------|
| FP        | CDK fingerprint                                                                    |
| ExtFP     | Extended CDK fingerprint (includes 25 bits for ring features and isotopic masses) |
| EStateFP  | Electrotopological state fingerprint (79 bits)                                    |
| GraphFP   | CDK fingerprinter ignoring bond orders                                            |
| MACCSFP   | Public MACCS fingerprint                                                          |
| PubchemFP | PubChem substructure fingerprint                                                  |
| SubFP     | Fingerprint describing 307 substructures                                          |
| KRFP      | Klekota-Roth fingerprint                                                          |
| AP2DFP    | Atom pair 2D fingerprint as implemented in PaDEL                                  |
| HybridFP  | CDK fingerprint ignoring aromaticity                                              |
| LingoFP   | LINGO fingerprint                                                                 |
| SPFP      | Fingerprint based on the shortest paths between two atoms                        |
| SigFP     | Signature fingerprint                                                             |
| CircFP    | Circular fingerprint                                                              |

### ⚡ Parallel processing

Speed things up by spreading molecules across several CPU cores. Each worker runs its own single-core-pinned JVM, so parallelism comes purely from the number of processes spawned — not from oversubscribing the host:

```python
cdk = CDK()
print(cdk.calculate(mols, njobs=8))
```

By default, molecules are auto-balanced evenly across `njobs` workers (`chunksize=None`), which minimizes JVM startup overhead while keeping every worker busy — the fastest setting for most workloads. A fixed `chunksize` can be provided instead if finer control is needed.

### 🔍 Details about descriptors

Details about all descriptors and fingerprints — name, description, type and dimensionality — can be obtained as a DataFrame:

```python
print(CDK.get_details())
```

Or for a single descriptor by name:

```python
print(CDK.get_details('MW'))
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/OlivierBeq/CDK_pywrapper/blob/master/LICENSE) file for details.

## 📚 API Documentation

```python
class CDK(ignore_3D=True, fingerprint=None, nbits=1024, depth=6, backend_smiles=False):
```

Constructor of a CDK calculator for molecular descriptors or a fingerprint.

#### Parameters

- ***ignore_3D  : bool***
  Whether to exclude the 65 3D molecular descriptors (default: `True`). Ignored if `fingerprint` is set.
- ***fingerprint  : FPType | None***
  Type of fingerprint to calculate (default: `None`). If `None`, descriptors are calculated instead.
- ***nbits  : int***
  Number of bits in the fingerprint, for fingerprints with a configurable size (default: 1024).
- ***depth  : int***
  Search depth of the fingerprint, for fingerprints with a configurable depth (default: 6).
- ***backend_smiles  : bool***
  Use SMILES rather than the V2000 SD format as the interchange format with the CDKdesc backend
  (default: `False`). Ignored if `ignore_3D=False`.

________________

```python
def calculate(mols, show_banner=True, cdk_smiles=False, njobs=1, chunksize=None):
```

Calculates CDK molecular descriptors and/or a fingerprint. Installs a matching JRE on first use if none is found.

#### Parameters

- ***mols  : Iterable[Chem.Mol]***
  RDKit molecule objects for which to obtain CDK descriptors/fingerprints (must have 3D
  conformers if 3D descriptors are requested).
- ***show_banner  : bool***
  Displays default notice about CDK.
- ***cdk_smiles  : bool***
  If `True`, also return the canonical SMILES CDK generated for each parsed molecule.
- ***njobs  : int***
  Number of concurrent processes used to calculate descriptors/fingerprints in parallel; must
  not exceed the number of available CPU cores. Each spawned Java process is pinned to a single
  core (`-XX:ActiveProcessorCount=1`), since parallelism comes from spawning `njobs` OS processes
  rather than from letting each JVM oversubscribe the host's full core count.
- ***chunksize  : int | None***
  Number of molecules processed per worker process. If `None` (default), molecules are
  auto-balanced across all `njobs` workers so every worker gets work. Ignored if `njobs` is 1.
- ***return_type  : pd.DataFrame***
  Pandas DataFrame containing CDK molecular descriptors and/or fingerprint values, one row per molecule.

________________

```python
@staticmethod
def get_details(desc_name=None):
```

Obtain metadata about either one or all descriptors/fingerprint bits.

#### Parameters

- ***desc_name  : str | None***
  Name of the descriptor to obtain details about (default: `None`). If `None`, returns details
  about all 288 descriptors.
- ***return_type  : pd.DataFrame***
  Pandas DataFrame with columns `Name`, `Description`, `Type` and `Dimensions`.
