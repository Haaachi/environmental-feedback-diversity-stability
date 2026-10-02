"""Regression checks for manuscript PR across both model storage formats."""

import importlib.util
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/plot/plot_main_stability_diversity_phase.py"
SPEC = importlib.util.spec_from_file_location("phase_maps", SCRIPT)
PHASE_MAPS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PHASE_MAPS)


class ManuscriptPRTests(unittest.TestCase):
    def compute(self, b10, legacy=False):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            if legacy:
                with h5py.File(root / "CR_phase2_meta.h5", "w") as h5:
                    h5["rmean_list"] = [1.0]
                    h5["hstress_list"] = [0.1]
                    h5.attrs.update(n_ax1=1, n_ax2=1, n_ensemble=b10.shape[1], S=b10.shape[2])
                rows = root / "phase2_rows"
                rows.mkdir()
                with h5py.File(rows / "phase2_row_001.h5", "w") as h5:
                    h5["B_last3"] = b10[..., -3:]
                compute = PHASE_MAPS.compute_maps_legacy_phase2
            else:
                with h5py.File(root / "phase_r_h_meta.h5", "w") as h5:
                    h5["axis1"] = [1.0]
                    h5["axis2"] = [0.1]
                    h5.attrs["S"] = b10.shape[2]
                rows = root / "raw/phase_r_h_rows"
                rows.mkdir(parents=True)
                with h5py.File(rows / "phase_r_h_row_0001.h5", "w") as h5:
                    h5["B_last10"] = b10
                    h5["completed_axis2"] = [True]
                compute = PHASE_MAPS.compute_maps
            return compute(root, b_thresh=1e-3, cv_thresh=0.1, collapse_thresh=1e-3, min_group=10)

    def test_earlier_passages_do_not_change_manuscript_pr(self):
        b10 = np.array([
            [2., 4, 2, 4, 2, 4, 2, 2, 4, 2],
            [4., 2, 4, 2, 4, 2, 4, 2, 2, 2],
        ])[None, None, ...]
        original = self.compute(b10)
        changed = b10.copy()
        changed[..., :7] *= np.array([1., 10.])[None, None, :, None]
        perturbed = self.compute(changed)
        for key in (
            "std_participation_ratio_B",
            "std_participation_ratio_B_osc_mean",
            "std_participation_ratio_B_osc_norm_mean",
        ):
            np.testing.assert_allclose(original[key], perturbed[key])
        np.testing.assert_allclose(original["std_participation_ratio_B"], 1.0)
        np.testing.assert_allclose(original["std_participation_ratio_B_osc_norm_mean"], 0.5)

    def test_small_species_variation_is_retained_in_both_formats(self):
        last3 = np.array([[0.01, 0.02, 0.03], [0.0001, 0.0002, 0.0003]])
        b10 = np.concatenate([np.repeat(last3[:, :1], 7, axis=-1), last3], axis=-1)[None, None, ...]
        expected_pr = (1.0 + 0.01) ** 2 / (1.0 + 0.01 ** 2)
        current = self.compute(b10)
        legacy = self.compute(b10, legacy=True)
        for result in (current, legacy):
            np.testing.assert_allclose(result["std_participation_ratio_B"], expected_pr, rtol=1e-6)
            np.testing.assert_allclose(result["std_participation_ratio_B_osc_norm_mean"], expected_pr / 2, rtol=1e-6)
        for key in ("community_CV_B", "is_oscillating_B", "is_collapsed_B", "std_participation_ratio_B"):
            np.testing.assert_allclose(current[key], legacy[key])

    def test_no_temporal_variation_has_undefined_participation(self):
        result = PHASE_MAPS.std_participation_ratio(np.ones((2, 3)))
        self.assertTrue(np.isnan(result))


if __name__ == "__main__":
    unittest.main()
