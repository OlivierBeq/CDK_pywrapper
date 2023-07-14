# -*- coding: utf-8 -*-


"""Tests for molecular fingerprints."""

import unittest

from CDK_pywrapper import CDK, FPType
from tests.constants import *


class TestFingerprints(unittest.TestCase):
    """Tests for CDK_pywrapper molecular fingerprints."""

    def setUp(self) -> None:
        self.molecules = list(MOLECULES.values())

    def test_fingerprint(self):
        sizes = {FPType.EStateFP: 79, FPType.MACCSFP: 166, FPType.PubchemFP: 881,
                 FPType.KRFP: 4860, FPType.SubFP: 307, FPType.AP2DFP: 780}
        for fp_type in FPType:
            cdk = CDK(fingerprint=fp_type)
            values = cdk.calculate(self.molecules, show_banner=False)
            if fp_type is not FPType.SigFP:
                self.assertEqual(values.shape, (len(MOLECULES), sizes.get(fp_type, 1024)))
                self.assertEqual(len(values.columns.unique().tolist()), sizes.get(fp_type, 1024))
            else:
                self.assertTrue(values.shape[0] == len(MOLECULES))
                self.assertEqual(len(values.columns.unique().tolist()), len(values.columns))
            self.assertFalse(values.isna().any().any())


    def test_fingerprint_multithread(self):
        sizes = {FPType.EStateFP: 79, FPType.MACCSFP: 166, FPType.PubchemFP: 881,
                 FPType.KRFP: 4860, FPType.SubFP: 307, FPType.AP2DFP: 780}
        for fp_type in FPType:
            cdk = CDK(fingerprint=fp_type)
            values = cdk.calculate(self.molecules, show_banner=False, njobs=-1, chunksize=1)
            if fp_type is not FPType.SigFP:
                self.assertEqual(values.shape, (len(MOLECULES), sizes.get(fp_type, 1024)))
                self.assertEqual(len(values.columns.unique().tolist()), sizes.get(fp_type, 1024))
            else:
                self.assertTrue(values.shape[0] == len(MOLECULES))
                self.assertEqual(len(values.columns.unique().tolist()), len(values.columns))
            self.assertFalse(values.isna().any().any())
