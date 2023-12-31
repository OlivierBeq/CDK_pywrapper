# -*- coding: utf-8

"""Python wrapper for CDK descriptors and fingerprints"""

from __future__ import annotations

import io
import multiprocessing
import os
import subprocess
import warnings
from copy import deepcopy
from enum import Enum, auto
from subprocess import Popen, PIPE
from typing import Iterable, List, Optional

import more_itertools
import numpy as np
import pandas as pd
from bounded_pool_executor import BoundedProcessPoolExecutor
from rdkit import Chem
from rdkit.rdBase import BlockLogs

from .utils import install_java, mktempfile, needsHs


class FPType(Enum):
    FP = auto()
    ExtFP = auto()
    EStateFP = auto()
    GraphFP = auto()
    MACCSFP = auto()
    PubchemFP = auto()
    SubFP = auto()
    KRFP = auto()
    AP2DFP = auto()
    HybridFP = auto()
    LingoFP = auto()
    SPFP = auto()
    SigFP = auto()
    CircFP = auto()


class CDK:
    """Wrapper to obtain molecular descriptor from CDK."""

    lock = multiprocessing.RLock()  # Ensure installation of JRE is thread safe
    # Path to the JAR file
    _jarfile = os.path.abspath(os.path.join(__file__, os.pardir, 'CDKdesc.jar'))

    def __init__(self, ignore_3D: bool = True, fingerprint: FPType = None, nbits: int = 1024, depth: int = 6):
        """Instantiate a wrapper to calculate CDK molecular descriptors or a fingerprint.

        :param ignore_3D: whether to include 3D molecular descriptors
        :param fingerprint: a fingerprint type to be calculated (default: None, calculates descriptors)
        :param nbits: number of bits (default: 1024 unless the fingerprint has a fixed size)
        :param depth: depth of the fingerprint (default: 6 unless the fingerprint does not depend on depth)
        """
        # Ensure the jar file exists
        if not os.path.isfile(self._jarfile):
            raise IOError('The required CDKdesc JAR file is not present. Reinstall CDK-pywrapper.')
        if fingerprint is not None:
            if not isinstance(fingerprint, FPType):
                raise TypeError(f'Fingerprint type not supported: {fingerprint}')
        self.include_3D = not ignore_3D
        self.fingerprint = None if fingerprint is None else fingerprint.name
        self.nbits = nbits
        self.depth = depth

    def calculate(self, mols: List[Chem.Mol], show_banner: bool = True, njobs: int = 1,
                  chunksize: Optional[int] = 1000) -> pd.DataFrame:
        """Calculate molecular fingerprints.

        :param mols: RDKit molecules for which descriptors/fingerprints should be calculated
        (must have 3D conformers if calculating descriptors)
        :param show_banner: If True, show notice on this package usage
        :param njobs: number of concurrent processes
        :param chunksize: number of molecules to be processed by a process; ignored if njobs is 1
        :return: a pandas DataFrame containing all CDK descriptor/fingerprint values
        """
        if show_banner:
            self._show_banner()
        if njobs < 0:
            njobs = os.cpu_count() - njobs + 1
        # Parallelize should need be
        if njobs > 1:
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                futures = [worker.submit(self._multiproc_calculate, list(chunk))
                           for chunk in more_itertools.batched(mols, chunksize)
                           ]
            return (pd.concat([future.result() for future in futures]).
                    reset_index(drop=True)
                    .convert_dtypes()
                    )
        # Single process
        return self._calculate(list(mols))

    def _show_banner(self):
        """Print info message for citing."""
        print("""The Chemistry Development Kit (CDK) is a collection of modular Java libraries
for processing chemical information (Cheminformatics). It can compute 14 different fingerprint
types and 287 molecular descriptors (it requires 3D molecular structures for the latter).

###################################

Should you publish results based on the PaDEL descriptors,
please cite:

Willighagen et al., (2017)  J. Cheminf. 9(3), doi:10.1186/s13321-017-0220-4,
May and Steinbeck., (2014) J. Cheminf., doi:10.1186/1758-2946-6-3,
Steinbeck et al., (2006) Curr. Pharm. Des. 12(17):2111-2120, doi:10.2174/138161206777585274,
Steinbeck et al., (2003) J. Chem. Inf. Comput. Sci. 43(2):493-500, doi:10.1021/ci025584y.

###################################

""")

    def _prepare_command(self, mols: List[Chem.Mol]) -> str:
        """Create the CDK command to be run to obtain molecular descriptors.

        :param mols: molecules to obtained molecular descriptors of
        :return: The command to run.
        """
        # 1) Ensure JRE is accessible
        with self.lock:
            self._java_path = install_java(19)
        # 2) Create temp SD v2k file
        self._tmp_sd = mktempfile('molecules_v2k.sd')
        self._n_mols = 0
        self._skipped = []
        self.n = 0
        try:
            block = BlockLogs()
            writer = Chem.SDWriter(self._tmp_sd)
            # Ensure V2000 as CDK cannot properly process v3000
            writer.SetForceV3000(False)
            for i, mol in enumerate(mols):
                if mol is not None and isinstance(mol, Chem.Mol):
                    if mol.GetNumAtoms() > 999:
                        raise ValueError('Cannot calculate descriptors for molecules with more than 999 atoms.')
                    # Does molecule lack hydrogen atoms?
                    if needsHs(mol):
                        warnings.warn('Molecule lacks hydrogen atoms: this might affect the value of calculated descriptors')
                    # Are molecules 3D
                    if self.include_3D:
                        confs = list(mol.GetConformers())
                        if self.fingerprint is None and not (len(confs) > 0 and confs[-1].Is3D()):
                            raise ValueError('Cannot calculate the 3D descriptors of a conformer-less molecule')
                    writer.write(mol)
                    self._n_mols += 1
                else:
                    self._skipped.append(i)
                self.n += 1
            writer.close()
            del block
        except ValueError as e:
            # Free resources and raise error
            writer.close()
            del block
            os.remove(self._tmp_sd)
            raise e from None
        # 3) Create command
        java_path = install_java(19)
        command_parameters = (f"-f {self.fingerprint} -nBits {self.nbits} "
                              f"-depth {self.depth}") if self.fingerprint is not None else ""
        command = f"{java_path} -jar {self._jarfile} -i {self._tmp_sd} {command_parameters}"
        return command

    def _cleanup(self) -> None:
        """Cleanup resources used for calculation."""
        # Remove temporary files
        os.remove(self._tmp_sd)

    def _run_command(self, command: str) -> pd.DataFrame:
        """Run the CDK command.

        :param command: The command to be run.
        """
        with Popen(command.split(), stdout=PIPE, stderr=subprocess.DEVNULL) as process:
            values = process.stdout.read().decode()
        # CDK barf preventing correct parsing
        if 'not found' in values:
            # Omit error
            values = '\n'.join(line for line in values.split('\n') if 'not found' not in line)
        # Empty result file
        if len(values) == 0:
            details = self.get_details()
            values = pd.DataFrame(np.full((self._n_mols, details.shape[0]), np.nan),
                                  columns=details.Name)
        elif '{' not in values:
            values = pd.read_csv(io.StringIO(values), sep=' ')
        else:
            try:
                values = pd.DataFrame.from_dict(eval('{%s}' % values), orient='index').fillna(0)
            except pd.errors.EmptyDataError:
                raise RuntimeError('CDK could not obtain molecular descriptors, maybe due to a faulty molecule')
        # If only 2D, remove 3D descriptors
        if not self.include_3D and self.fingerprint is None:
            # Get 3D descriptor names to remove
            descs_3D = self.get_details()
            descs_3D = descs_3D[descs_3D.Dimensions == '3D']
            values = values[[col for col in values.columns if col not in descs_3D.Name.tolist()]]
        return values

    def _calculate(self, mols: List[Chem.Mol]) -> pd.DataFrame:
        """Calculate CDK molecular descriptors on one process.

        :param mols: RDKit molecules for which CDK descriptors and fingerprints should be calculated.
        :return: a pandas DataFrame containing CDK descriptor values
        """
        # Prepare inputs
        command = self._prepare_command(mols)
        # Run command and obtain results
        results = self._run_command(command)
        # Cleanup
        self._cleanup()
        # Insert lines of skipped molecules
        if len(self._skipped):
            results = pd.DataFrame(np.insert(results.values, self._skipped,
                                              values=[np.NaN] * len(results.columns),
                                              axis=0),
                                    columns=results.columns)
        results = (results.apply(pd.to_numeric, errors='coerce', axis=1)
                          .convert_dtypes()
                   )
        return results

    def _multiproc_calculate(self, mols: List[Chem.Mol]) -> pd.DataFrame:
        """Calculate CDK descriptors and fingerprints in thread-safe manner.

        :param mols: RDKit molecules for which CDK descriptors and fingerprints should be calculated
        :return: a pandas DataFrame containing all CDK descriptor values
        """
        # Copy self instance to make thread safe
        cdk = deepcopy(self)
        # Run copy
        result = cdk.calculate(mols, show_banner=False, njobs=1)
        return result

    @staticmethod
    def get_details(desc_name: Optional[str] = None):
        """Obtain details about either one or all descriptors.

        :param desc_name: the name of the descriptor to obtain details about (default: None).
        If None, returns details about all descriptors.
        """
        details = pd.read_json(os.path.abspath(os.path.join(__file__, os.pardir, 'descs.json')), orient='index')
        if desc_name is not None:
            if desc_name not in details.Name.tolist():
                raise ValueError(f'descriptor name {desc_name} is not available')
            details = details[details.Name == desc_name]
        return details
