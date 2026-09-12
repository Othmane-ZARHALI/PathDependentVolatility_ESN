"""Tests for the paired, checkpointed M2/M3/M2+O/M4 experiment."""

from __future__ import annotations

import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from esn_eq import (
    AcceptanceCriteria,
    ArchitectureConfig,
    FactorialDesign,
    FactorialExperimentConfig,
    FactorialExperimentRecorder,
    MarketEnvironment,
    ModelVariant,
    MultiSeedFactorialExperiment,
    SimulationConfig,
    StructuralParameters,
)
from esn_eq.config import BankConfig, EchoStateConfig
from esn_eq.features import FeatureEngine
from esn_eq.sampling import ParameterSpace
from esn_eq.reporting import FactorialMarkdownReport
from esn_eq.states import build_mode_layout


def factorial_architecture() -> ArchitectureConfig:
    """Return a small maximal architecture suitable for integration tests."""
    return ArchitectureConfig(
        feedback=BankConfig((1.0, 10.0), 0.80, True),
        spike=BankConfig((20.0,), 0.50, True),
        cluster=BankConfig((0.2, 2.0), 0.0, False),
        orthogonal=BankConfig((0.5, 5.0), 0.0, False),
        echo=EchoStateConfig(enabled=True, units=4, seed=101),
    )


def factorial_simulation() -> SimulationConfig:
    """Provide enough observations for every configured diagnostic window."""
    return SimulationConfig(
        years=3.0,
        burn_years=0.5,
        paths=3,
        observations_per_year=40,
        steps_per_observation=1,
    )


class FactorialDesignTests(unittest.TestCase):
    """Verify the exact-zero controls and architectural separation."""

    def test_cells_share_core_parameters_and_apply_only_optional_switches(self) -> None:
        """Every four-cell block must be a genuinely paired parameter draw."""
        design = FactorialDesign.sample(3, seed=202, method="sobol")
        for candidate_id in range(3):
            cells = [item for item in design.cells() if item.candidate_id == candidate_id]
            self.assertEqual([item.variant for item in cells], list(ModelVariant))
            core = {
                (
                    item.parameters.cluster_loading,
                    item.parameters.feedback_curvature,
                    item.parameters.feedback_shift,
                    item.parameters.spike_curvature,
                    item.parameters.spike_shift,
                    item.parameters.cap_level,
                    item.parameters.cap_sharpness,
                )
                for item in cells
            }
            self.assertEqual(len(core), 1)
            self.assertEqual(cells[0].parameters.echo_loading, 0.0)
            self.assertEqual(cells[0].parameters.orthogonal_curvature, 0.0)
            self.assertEqual(cells[1].parameters.orthogonal_curvature, 0.0)
            self.assertEqual(cells[2].parameters.echo_loading, 0.0)

    def test_custom_parameter_space_is_retained_for_the_manifest(self) -> None:
        """A targeted design must carry its sampling measure into audit records."""
        space = ParameterSpace.targeted_stage2()
        design = FactorialDesign.sample(
            2,
            seed=203,
            method="sobol",
            space=space,
            space_name="targeted_stage2",
        )
        self.assertEqual(design.parameter_space, space)
        self.assertEqual(design.space_name, "targeted_stage2")

    def test_standardized_orthogonal_energy_has_stationary_unit_scale(self) -> None:
        """Centered energy should be comparable in scale with other score features."""
        architecture = replace(
            factorial_architecture(),
            orthogonal_energy_mode="standardized_centered",
        )
        layout = build_mode_layout(architecture)
        generator = np.random.default_rng(305)
        states = generator.standard_normal((50000, 1, architecture.factor_count))
        energy = FeatureEngine(architecture, layout).build(states).orthogonal_energy
        self.assertAlmostEqual(float(np.mean(energy)), 0.0, delta=0.03)
        self.assertAlmostEqual(float(np.std(energy)), 1.0, delta=0.04)

    def test_orthogonal_bank_does_not_redefine_the_echo_feature(self) -> None:
        """M3 and M4 must use the same echo map on identical common-core states."""
        core = replace(ArchitectureConfig(), echo=EchoStateConfig(enabled=True, seed=303))
        maximal = replace(core, orthogonal=BankConfig((0.5, 2.0)))
        core_layout = build_mode_layout(core)
        maximal_layout = build_mode_layout(maximal)
        generator = np.random.default_rng(304)
        core_states = generator.standard_normal((12, 3, core.factor_count))
        maximal_states = generator.standard_normal((12, 3, maximal.factor_count))
        maximal_states[..., : core.factor_count] = core_states
        core_echo = FeatureEngine(core, core_layout).build(core_states).echo
        maximal_echo = FeatureEngine(maximal, maximal_layout).build(maximal_states).echo
        np.testing.assert_array_equal(core_echo, maximal_echo)


class FactorialExperimentTests(unittest.TestCase):
    """Exercise multi-seed pairing, summaries, checkpoints, and reproducibility."""

    def test_checkpointed_run_is_complete_and_reproducible(self) -> None:
        """Two cache seeds must yield four matched cells for every core draw."""
        runner = MultiSeedFactorialExperiment(
            factorial_architecture(),
            factorial_simulation(),
            MarketEnvironment(),
            AcceptanceCriteria.illustrative(),
        )
        design = FactorialDesign.sample(2, seed=401, method="sobol")
        config = FactorialExperimentConfig((402, 403))
        structural = (StructuralParameters(),)
        with tempfile.TemporaryDirectory() as directory:
            recorder = FactorialExperimentRecorder(directory)
            first = runner.run(structural, design, config, recorder)
            second = runner.run(structural, design, config)
            pd.testing.assert_frame_equal(first.runs, second.runs)
            self.assertEqual(len(first.runs), 16)
            self.assertEqual(len(first.seed_summary), 8)
            self.assertFalse(first.scenario_summary.empty)
            self.assertFalse(first.criterion_summary.empty)
            self.assertTrue((Path(directory) / "checkpoints" / "seed_402.csv").is_file())
            self.assertTrue((Path(directory) / "checkpoints" / "seed_403.csv").is_file())
            self.assertTrue((Path(directory) / "runs.csv").is_file())
            report = FactorialMarkdownReport().build(first)
            self.assertIn("Acceptance by independent cache seed", report)
            self.assertIn("Individual gate pass rates", report)
            self.assertIn("Post-hoc threshold sensitivity", report)
            self.assertFalse(FactorialMarkdownReport.stress_summary(first.runs).empty)
        counts = first.runs.groupby(["seed", "candidate_id"])["variant"].nunique()
        self.assertTrue((counts == 4).all())
        self.assertIn("echo_orthogonal_interaction", set(first.paired_effects.effect))

    def test_factorial_runner_rejects_echo_dependence_on_orthogonal_states(self) -> None:
        """An entangled echo would invalidate the intended two-by-two contrast."""
        architecture = factorial_architecture()
        entangled = replace(
            architecture,
            echo=replace(
                architecture.echo,
                input_banks=("feedback", "spike", "cluster", "orthogonal"),
            ),
        )
        with self.assertRaises(ValueError):
            MultiSeedFactorialExperiment(
                entangled,
                factorial_simulation(),
                MarketEnvironment(),
                AcceptanceCriteria.illustrative(),
            )


if __name__ == "__main__":
    unittest.main()
