"""Tests for parallel dispatch."""

import unittest

from CDK_pywrapper import CDK
from CDK_pywrapper.cdk_pywrapper import _make_chunks
from tests.constants import MOLECULES, NJOBS


class TestParallel(unittest.TestCase):
    """Tests for CDK_pywrapper parallel dispatch."""

    def setUp(self) -> None:
        self.molecules = list(MOLECULES.values())

    def test_make_chunks_balances_when_chunksize_is_none(self):
        chunks = _make_chunks(self.molecules, njobs=3, chunksize=None)
        self.assertEqual(len(chunks), min(3, len(self.molecules)))
        sizes = sorted(len(chunk) for chunk in chunks)
        self.assertLessEqual(sizes[-1] - sizes[0], 1)
        self.assertEqual(sum(sizes), len(self.molecules))

    def test_make_chunks_never_exceeds_number_of_molecules(self):
        chunks = _make_chunks(self.molecules, njobs=100, chunksize=None)
        self.assertEqual(len(chunks), len(self.molecules))

    def test_make_chunks_respects_explicit_chunksize(self):
        chunks = _make_chunks(self.molecules, njobs=2, chunksize=2)
        self.assertTrue(all(len(chunk) <= 2 for chunk in chunks))

    def test_make_chunks_empty_input(self):
        self.assertEqual(_make_chunks([], njobs=4, chunksize=None), [])

    def test_njobs_below_one_raises_valueerror(self):
        cdk = CDK()
        with self.assertRaises(ValueError):
            cdk.calculate(self.molecules, show_banner=False, njobs=0)

    def test_njobs_above_cpu_count_raises_valueerror(self):
        import os
        cdk = CDK()
        with self.assertRaises(ValueError):
            cdk.calculate(self.molecules, show_banner=False, njobs=(os.cpu_count() or 1) + 1)

    def test_parallel_matches_single_process(self):
        cdk = CDK()
        single = cdk.calculate(self.molecules, show_banner=False, njobs=1)
        parallel = cdk.calculate(self.molecules, show_banner=False, njobs=NJOBS, chunksize=None)
        self.assertEqual(single.shape, parallel.shape)
        self.assertListEqual(single.columns.tolist(), parallel.columns.tolist())


if __name__ == '__main__':
    unittest.main()
