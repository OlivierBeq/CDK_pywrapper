# -*- coding: utf-8 -*-

"""Tests for molecular descriptors."""

import unittest

import numpy as np

from CDK_pywrapper import CDK
from tests.constants import *


class TestDescriptors(unittest.TestCase):
    """Tests for CDK_pywrapper molecular descriptors."""

    def setUp(self) -> None:
        """Create the molecular descriptor calculator."""
        self.cdk = CDK()
        self.cdk3d = CDK(ignore_3D=False)
        self.molecules = list(MOLECULES.values())

    def test_2D_descriptor_size(self):
        values = self.cdk.calculate(self.molecules, show_banner=False)
        self.assertEqual(values.shape, (len(MOLECULES), 223))
        self.assertEqual(len(values.columns.unique().tolist()), 223)
        # Include SMILES
        values = self.cdk.calculate(self.molecules, show_banner=False, cdk_smiles=True)
        self.assertEqual(values.shape, (len(MOLECULES), 224))
        self.assertEqual(len(values.columns.unique().tolist()), 224)
        self.assertTrue(values.columns[0] == 'SMILES')
        self.assertTrue(values.iloc[:, 0].dtype == np.dtype('O'))

    def test_2D_descriptor_multithread(self):
        values = self.cdk.calculate(self.molecules, show_banner=False, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES), 223))
        self.assertEqual(len(values.columns.unique().tolist()), 223)
        # Include SMILES
        values = self.cdk.calculate(self.molecules, show_banner=False, cdk_smiles=True, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES), 224))
        self.assertEqual(len(values.columns.unique().tolist()), 224)
        self.assertTrue(values.columns[0] == 'SMILES')
        self.assertTrue(values.iloc[:, 0].dtype == np.dtype('O'))

    def test_2D_cisplatin(self):
        values = self.cdk.calculate([MOLECULES['cisplatin']], show_banner=False)
        self.assertEqual(values.shape, (1, 223))
        self.assertEqual(len(values.columns.unique().tolist()), 223)
        # Include SMILES
        values = self.cdk.calculate([MOLECULES['cisplatin']], show_banner=False, cdk_smiles=True)
        self.assertEqual(values.shape, (1, 224))
        self.assertEqual(len(values.columns.unique().tolist()), 224)
        self.assertTrue(values.columns[0] == 'SMILES')
        self.assertTrue(values.iloc[:, 0].dtype == np.dtype('O'))

    def test_3D_descriptor_size(self):
        values = self.cdk3d.calculate(self.molecules, show_banner=False)
        self.assertEqual(values.shape, (len(MOLECULES), 288))
        self.assertEqual(len(values.columns.unique().tolist()), 288)
        # Include SMILES
        values = self.cdk3d.calculate(self.molecules, show_banner=False, cdk_smiles=True)
        self.assertEqual(values.shape, (len(MOLECULES), 289))
        self.assertEqual(len(values.columns.unique().tolist()), 289)
        self.assertTrue(values.columns[0] == 'SMILES')
        self.assertTrue(values.iloc[:, 0].dtype == np.dtype('O'))

    def test_3D_descriptor_multithread(self):
        values = self.cdk3d.calculate(self.molecules, show_banner=False, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES), 288))
        self.assertEqual(len(values.columns.unique().tolist()), 288)
        # Include SMILES
        values = self.cdk3d.calculate(self.molecules, show_banner=False, cdk_smiles=True, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES), 289))
        self.assertEqual(len(values.columns.unique().tolist()), 289)
        self.assertTrue(values.columns[0] == 'SMILES')
        self.assertTrue(values.iloc[:, 0].dtype == np.dtype('O'))

    def test_3D_cisplatin(self):
        values = self.cdk3d.calculate([MOLECULES['cisplatin']], show_banner=False)
        self.assertEqual(values.shape, (1, 288))
        self.assertEqual(len(values.columns.unique().tolist()), 288)
        # Include SMILES
        values = self.cdk3d.calculate([MOLECULES['cisplatin']], show_banner=False, cdk_smiles=True)
        self.assertEqual(values.shape, (1, 289))
        self.assertEqual(len(values.columns.unique().tolist()), 289)
        self.assertTrue(values.columns[0] == 'SMILES')
        self.assertIsInstance(values.iloc[0, 0], str)

    def test_get_details(self):
        details = self.cdk.get_details()
        self.assertEqual(details.shape, (288, 4))
        self.assertListEqual(details.columns.tolist(), ['Name', 'Description', 'Type', 'Dimensions'])

    def test_2D_descriptor_smiles_backend(self):
        cdk = CDK(backend_smiles=True)
        values = cdk.calculate(self.molecules, show_banner=False)
        self.assertEqual(values.shape, (len(MOLECULES), 223))
        self.assertEqual(len(values.columns.unique().tolist()), 223)
        # Include SMILES
        values = cdk.calculate(self.molecules, show_banner=False, cdk_smiles=True)
        self.assertEqual(values.shape, (len(MOLECULES), 224))
        self.assertEqual(len(values.columns.unique().tolist()), 224)
        self.assertTrue(values.columns[0] == 'SMILES')
        self.assertTrue(values.iloc[:, 0].dtype == np.dtype('O'))

    def test_2D_descriptor_smiles_backend_multithread(self):
        cdk = CDK(backend_smiles=True)
        values = cdk.calculate(self.molecules, show_banner=False, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES), 223))
        self.assertEqual(len(values.columns.unique().tolist()), 223)
        # Include SMILES
        values = cdk.calculate(self.molecules, show_banner=False, cdk_smiles=True, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES), 224))
        self.assertEqual(len(values.columns.unique().tolist()), 224)
        self.assertTrue(values.columns[0] == 'SMILES')
        self.assertTrue(values.iloc[:, 0].dtype == np.dtype('O'))

    def test_3D_descriptor_smiles_backend(self):
        cdk = CDK(ignore_3D=False, backend_smiles=True)
        values = cdk.calculate(self.molecules, show_banner=False)
        self.assertEqual(values.shape, (len(MOLECULES), 288))
        self.assertEqual(len(values.columns.unique().tolist()), 288)
        # Include SMILES
        values = cdk.calculate(self.molecules, show_banner=False, cdk_smiles=True)
        self.assertEqual(values.shape, (len(MOLECULES), 289))
        self.assertEqual(len(values.columns.unique().tolist()), 289)
        self.assertTrue(values.columns[0] == 'SMILES')
        self.assertTrue(values.iloc[:, 0].dtype == np.dtype('O'))

    def test_3D_descriptor_smiles_backend_multithread(self):
        cdk = CDK(ignore_3D=False, backend_smiles=True)
        values = cdk.calculate(self.molecules, show_banner=False, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES), 288))
        self.assertEqual(len(values.columns.unique().tolist()), 288)
        # Include SMILES
        values = cdk.calculate(self.molecules, show_banner=False, cdk_smiles=True, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES), 289))
        self.assertEqual(len(values.columns.unique().tolist()), 289)
        self.assertTrue(values.columns[0] == 'SMILES')
        self.assertTrue(values.iloc[:, 0].dtype == np.dtype('O'))
