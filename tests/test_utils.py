"""Tests for utility functions."""

import os
import unittest
from unittest.mock import patch

from rdkit import Chem

from CDK_pywrapper.utils import get_java_in_dir, install_java, mktempdir, mktempfile, needsHs


class TestUtils(unittest.TestCase):
    """Tests for CDK_pywrapper utility functions."""

    def test_install_java(self):
        path = install_java(19)
        self.assertIsNotNone(path)
        self.assertTrue(os.path.isfile(path))

    def test_get_java_in_dir_missing_version_returns_none(self):
        self.assertIsNone(get_java_in_dir(os.path.dirname(__file__), 999999))

    def test_install_java_triggers_install_once_on_cold_cache(self):
        with patch('CDK_pywrapper.utils.get_java_in_dir', side_effect=[None, '/fake/java']) as mock_get, \
             patch('CDK_pywrapper.utils._jre_install') as mock_install:
            path = install_java(19)
        self.assertEqual(path, '/fake/java')
        mock_install.assert_called_once()
        self.assertEqual(mock_get.call_count, 2)

    def test_install_java_skips_install_on_warm_cache(self):
        with patch('CDK_pywrapper.utils.get_java_in_dir', return_value='/fake/java') as mock_get, \
             patch('CDK_pywrapper.utils._jre_install') as mock_install:
            path = install_java(19)
        self.assertEqual(path, '/fake/java')
        mock_install.assert_not_called()
        self.assertEqual(mock_get.call_count, 1)

    def test_mktempdir(self):
        path = mktempdir()
        self.assertTrue(os.path.isdir(path))

    def test_mktempfile(self):
        path = mktempfile()
        self.assertTrue(os.path.isfile(path))

    def test_needsHs_true_without_explicit_hydrogens(self):
        mol = Chem.MolFromSmiles('CCO')
        self.assertTrue(needsHs(mol))

    def test_needsHs_false_with_explicit_hydrogens(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
        self.assertFalse(needsHs(mol))


if __name__ == '__main__':
    unittest.main()
