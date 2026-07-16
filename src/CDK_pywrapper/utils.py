"""Utility functions."""

from __future__ import annotations

import glob
import os
import sys
import tempfile

from filelock import FileLock
from jdk import _JRE_DIR
from jdk import install as _jre_install
from rdkit import Chem


def install_java(version: int = 19) -> str | None:
    """Install a Java Runtime Environment.

    A cross-process file lock guards the check-then-install so that concurrent
    processes (e.g. parallel ``njobs`` workers, or pytest-xdist workers) racing on a
    cold JRE cache cannot collide on ``install-jdk``'s fixed, version-derived
    download/extraction path.

    :param version: the major JDK version to install (default: 19)
    :return: path to the ``java`` executable, or None if it could not be installed
    """
    os.makedirs(_JRE_DIR, exist_ok=True)
    with FileLock(os.path.join(_JRE_DIR, f'.install-{version}.lock')):
        path = get_java_in_dir(_JRE_DIR, version)
        if path is None:
            # Could not find JRE, install it
            _jre_install(version, jre=True)
            path = get_java_in_dir(_JRE_DIR, version)
    return path


def get_java_in_dir(dir: str, version: int) -> str | None:
    """Recursively search the directory to find a Java executable of the given version."""
    paths = glob.glob(os.path.join(dir, '**', 'bin',
                                  'java.exe' if sys.platform == "win32" else 'java'
                                  ), recursive=True)
    path = [path for path in paths if f'jdk-{version}' in path]
    if len(path):
        return os.path.abspath(path[0])
    return None


def mktempdir(suffix: str | None = None) -> str:
    """Return the path to a writeable temporary directory."""
    dir = tempfile.mkdtemp(suffix=suffix)
    return dir


def mktempfile(suffix: str | None = None) -> str:
    """Return the path to a writeable temporary file."""
    file = tempfile.mkstemp(suffix=suffix)
    os.close(file[0])
    return file[1]


def needsHs(mol: Chem.Mol) -> bool:
    """Return if the molecule lacks hydrogen atoms or not.

    :param mol: RDKit Molecule
    :return: True if the molecule lacks hydrogens.
    """
    for atom in mol.GetAtoms():
        nHNbrs = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                nHNbrs += 1
        noNeighbors = False
        if atom.GetTotalNumHs(noNeighbors) > nHNbrs:
            return True
    return False
