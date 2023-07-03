# -*- coding: utf-8 -*-
"""Python wrapper for CDK molecular descriptors and fingerprints."""

from .cdk_pywrapper import CDK, FPType  # noqa: F401 : imported but not used
from .version import get_version as __get_version__

__version__ = __get_version__()
