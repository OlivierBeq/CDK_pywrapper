<!--
<p align="center">
  <img src="https://github.com/OlivierBeq/CDK_pywrapper/raw/main/docs/source/logo.png" height="150">
</p>
-->

<h1 align="center">
  CDK_pywrapper
</h1>

<p align="center">
    <a href="https://github.com/OlivierBeq/CDK_pywrapper/actions/workflows/tests.yml">
        <img alt="Tests" src="https://github.com/OlivierBeq/CDK_pywrapper/workflows/Tests/badge.svg" />
    </a>
    <a href="https://pypi.org/project/CDK_pywrapper">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/CDK_pywrapper" />
    </a>
    <a href="https://pypi.org/project/CDK_pywrapper">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/CDK_pywrapper" />
    </a>
    <a href="https://github.com/OlivierBeq/CDK_pywrapper/blob/main/LICENSE">
        <img alt="PyPI - License" src="https://img.shields.io/pypi/l/CDK_pywrapper" />
    </a>
    <a href='https://CDK_pywrapper.readthedocs.io/en/latest/?badge=latest'>
        <img src='https://readthedocs.org/projects/CDK_pywrapper/badge/?version=latest' alt='Documentation Status' />
    </a>
    <a href="https://codecov.io/gh/OlivierBeq/CDK_pywrapper/branch/main">
        <img src="https://codecov.io/gh/OlivierBeq/CDK_pywrapper/branch/main/graph/badge.svg" alt="Codecov status" />
    </a>  
    <a href="https://github.com/cthoyt/cookiecutter-python-package">
        <img alt="Cookiecutter template from @cthoyt" src="https://img.shields.io/badge/Cookiecutter-snekpack-blue" /> 
    </a>
    <a href='https://github.com/psf/black'>
        <img src='https://img.shields.io/badge/code%20style-black-000000.svg' alt='Code style: black' />
    </a>
    <a href="https://github.com/OlivierBeq/CDK_pywrapper/blob/main/.github/CODE_OF_CONDUCT.md">
        <img src="https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg" alt="Contributor Covenant"/>
    </a>
</p>

Python wrapper for CDK molecular descriptors and fingerprints

## 💪 Getting Started

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
# Ensure hydrogens are explicit
mols = [Chem.AddHs(Chem.MolFromSmiles(smiles)) for smiles in smiles_list]
# Ensure molecules have 3D conformations
for mol in mols:
    _ = AllChem.EmbedMolecule(mol)

cdk = CDK()
print(cdk.calculate(mols))
```

The above calculates 287 molecular descriptors (23 1D, 200 2D and 64 3D).<br/>
:warning: Molecules are required to have conformers for descriptors to be calculated.<br/>

To obtain molecular fingerprint, one can use the following:

```python
from CDK_pywrapper import CDK, FPType

cdk = CDK(FPType.PubchemFP, nbits=2048)
print(cdk.calculate(mols))
```

The following fingerprints can be calculated:

| FPType    | Fingerprint name                                                                  |
|-----------|-----------------------------------------------------------------------------------|
| FP        | CDK fingerprint                                                                   |
| ExtFP     | Extended CDK fingerprint (includes 25 bits for ring features and isotopic masses) |
| GraphFP   | CDK fingerprinter ignoring bond orders                                            |
| MACCSFP   | Public MACCS fingerprint                                                          |
| PubchemFP | PubChem substructure fingerprint                                                  |
| SubFP     | Fingerprint describing 307 substructures                                          |
| KRFP      | Klekota-Roth fingerprint                                                          |
| AP2DFP    | Atom pair 2D fingerprint as implemented in PaDEL                                  |
| HybridFP  | CDK fingerprint ignoring aromaticity                                              |
| LingoFP   | LINGO fingerprint                                                                 |
| SPFP      | Fingerprint based on the shortest paths between two atoms                         |
| SigFP     | Signature fingerprint                                                             |
| CircFP    | Circular fingerprint                                                              |

## 🚀 Installation

The most recent release can be installed from
[PyPI](https://pypi.org/project/CDK_pywrapper/) with:

```shell
$ pip install CDK_pywrapper
```

The most recent code and data can be installed directly from GitHub with:

```bash
$ pip install git+https://github.com/OlivierBeq/CDK_pywrapper.git
```

## 👐 Contributing

Contributions, whether filing an issue, making a pull request, or forking, are appreciated. See
[CONTRIBUTING.md](https://github.com/OlivierBeq/CDK_pywrapper/blob/master/.github/CONTRIBUTING.md) for more information on getting involved.

## 👋 Attribution

### 📖 Citation

1. CDK version 2.8 ![CDK-badge.svg](https://zenodo.org/badge/DOI/10.5281/zenodo.7079512.svg) 
2. Willighagen et al. The Chemistry Development Kit (CDK) v2.0: atom typing, depiction, molecular formulas, and substructure searching. J. Cheminform. 2017; 9(3), doi:10.1186/s13321-017-0220-4
3. May and Steinbeck. Efficient ring perception for the Chemistry Development Kit. J. Cheminform. 2014, doi:10.1186/1758-2946-6-3
4. Steinbeck et al. Recent Developments of the Chemistry Development Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics. Curr. Pharm. Des. 2006; 12(17):2111-2120, doi:10.2174/138161206777585274 (free green Open Acccess version)
5. Steinbeck et al. The Chemistry Development Kit (CDK): An Open-Source Java Library for Chemo- and Bioinformatics. J. Chem. Inf. Comput. Sci. 2003 Mar-Apr; 43(2):493-500, doi:10.1021/ci025584y

### 🍪 Cookiecutter

This package was created with [@audreyfeldroy](https://github.com/audreyfeldroy)'s
[cookiecutter](https://github.com/cookiecutter/cookiecutter) package using [@cthoyt](https://github.com/cthoyt)'s
[cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack) template.

## 🛠️ For Developers

<details>
  <summary>See developer instructions</summary>

The final section of the README is for if you want to get involved by making a code contribution.

### Development Installation

To install in development mode, use the following:

```bash
$ git clone git+https://github.com/{{cookiecutter.github_organization_name}}/{{cookiecutter.github_repository_name}}.git
$ cd {{cookiecutter.github_repository_name}}
$ pip install -e .
```

### 🥼 Testing

After cloning the repository and installing `tox` with `pip install tox`, the unit tests in the `tests/` folder can be
run reproducibly with:

```shell
$ tox
```

Additionally, these tests are automatically re-run with each commit in a [GitHub Action](https://github.com/{{cookiecutter.github_organization_name}}/{{cookiecutter.github_repository_name}}/actions?query=workflow%3ATests).

### 📖 Building the Documentation

The documentation can be built locally using the following:

```shell
$ git clone git+https://github.com/{{cookiecutter.github_organization_name}}/{{cookiecutter.github_repository_name}}.git
$ cd {{cookiecutter.github_repository_name}}
$ tox -e docs
$ open docs/build/html/index.html
``` 

The documentation automatically installs the package as well as the `docs`
extra specified in the [`setup.cfg`](setup.cfg). `sphinx` plugins
like `texext` can be added there. Additionally, they need to be added to the
`extensions` list in [`docs/source/conf.py`](docs/source/conf.py).

### 📦 Making a Release

After installing the package in development mode and installing
`tox` with `pip install tox`, the commands for making a new release are contained within the `finish` environment
in `tox.ini`. Run the following from the shell:

```shell
$ tox -e finish
```

This script does the following:

1. Uses [Bump2Version](https://github.com/c4urself/bump2version) to switch the version number in the `setup.cfg`,
   `src/{{cookiecutter.package_name}}/version.py`, and [`docs/source/conf.py`](docs/source/conf.py) to not have the `-dev` suffix
2. Packages the code in both a tar archive and a wheel using [`build`](https://github.com/pypa/build)
3. Uploads to PyPI using [`twine`](https://github.com/pypa/twine). Be sure to have a `.pypirc` file configured to avoid the need for manual input at this
   step
4. Push to GitHub. You'll need to make a release going with the commit where the version was bumped.
5. Bump the version to the next patch. If you made big changes and want to bump the version by minor, you can
   use `tox -e bumpversion -- minor` after.
</details>