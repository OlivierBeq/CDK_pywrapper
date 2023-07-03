# -*- coding: utf-8 -*-
"""Tests for molecular descriptors."""

import unittest

from CDK_pywrapper import CDK
from tests.constants import MOLECULES


class TestDescriptors(unittest.TestCase):
    """Tests for CDK_pywrapper molecular descriptors."""
    def setUp(self) -> None:
        """Create the molecular descriptor calculator."""
        self.cdk = CDK()
        self.molecules = list(MOLECULES.values())

    def test_descriptor_size(self):
        """Test the dimensions of the output descriptor dataframe."""
        values = self.cdk.calculate(self.molecules, show_banner=False)
        self.assertEqual(values.shape, (len(MOLECULES), 287))
        self.assertFalse(values.isna().any().any())
        self.assertEqual(len(values.columns.unique().tolist()), 287)

    def test_descriptor_multithread(self):
        """Test the dimensions of the output descriptor dataframes calculated by different processes."""
        values = self.cdk.calculate(self.molecules, show_banner=False, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES), 287))
        self.assertFalse(values.isna().any().any())
        self.assertEqual(len(values.columns.unique().tolist()), 287)

    def test_cisplatin(self):
        """Test output of cisplatin."""
        values = self.cdk.calculate([MOLECULES["cisplatin"]], show_banner=False)
        self.assertNotEqual(values.shape, (len(MOLECULES), 287))
        self.assertFalse(values.isna().any().any())
        self.assertNotEqual(len(values.columns.unique().tolist()), 287)
