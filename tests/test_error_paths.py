"""Tests for error paths and edge cases."""

import unittest
from unittest.mock import patch

from rdkit import Chem

from CDK_pywrapper import CDK
from tests.constants import MOLECULES


class TestErrorPaths(unittest.TestCase):
    """Tests for CDK_pywrapper error paths."""

    def test_missing_jar_raises_oserror(self):
        with patch.object(CDK, '_jarfile', '/nonexistent/CDKdesc.jar'):
            with self.assertRaises(OSError):
                CDK()

    def test_invalid_fingerprint_type_raises_typeerror(self):
        with self.assertRaises(TypeError):
            CDK(fingerprint='not-a-fptype')

    def test_oversized_molecule_raises_valueerror(self):
        mol = Chem.MolFromSmiles('C' * 1005)
        self.assertGreater(mol.GetNumAtoms(), 999)
        cdk = CDK()
        with self.assertRaises(ValueError):
            cdk.calculate([mol], show_banner=False)

    def test_conformerless_molecule_raises_for_3d_descriptors(self):
        mol = Chem.MolFromSmiles('CCO')
        cdk = CDK(ignore_3D=False)
        with self.assertRaises(ValueError):
            cdk.calculate([mol], show_banner=False)

    def test_conformerless_molecule_succeeds_for_2d_descriptors(self):
        mol = Chem.MolFromSmiles('CCO')
        cdk = CDK()
        values = cdk.calculate([mol], show_banner=False)
        self.assertEqual(values.shape[0], 1)

    def test_missing_hydrogens_warns(self):
        mol = Chem.MolFromSmiles('CCO')
        cdk = CDK()
        with self.assertWarns(UserWarning):
            cdk.calculate([mol], show_banner=False)

    def test_none_molecule_skipped_and_reinserted_as_nan(self):
        cdk = CDK()
        mols = [MOLECULES['erlotinib'], None, MOLECULES['erlotinib']]
        values = cdk.calculate(mols, show_banner=False)
        self.assertEqual(values.shape[0], 3)
        self.assertTrue(values.iloc[1].isna().all())
        self.assertFalse(values.iloc[0].isna().any())
        self.assertFalse(values.iloc[2].isna().any())

    def test_get_details_invalid_name_raises_valueerror(self):
        with self.assertRaises(ValueError):
            CDK.get_details('this-descriptor-does-not-exist')

    def test_get_details_valid_name(self):
        details = CDK.get_details('MW')
        self.assertEqual(details.shape[0], 1)
        self.assertEqual(details.Name.iloc[0], 'MW')


if __name__ == '__main__':
    unittest.main()
