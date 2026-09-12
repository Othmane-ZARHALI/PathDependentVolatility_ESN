"""Human-readable reports for paired factorial experiments."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from .factorial import FactorialExperimentResult, VARIANT_ORDER


DISPLAY_METRICS = (
    ("max_abs_return_acf", "max |return ACF|"),
    ("mean_log_variance_acf", "mean log-variance ACF"),
    ("rough_hurst", "rough H"),
    ("leverage_rank", "rank leverage"),
    ("zumbach_rank", "rank Zumbach"),
    ("excess_kurtosis", "excess kurtosis"),
    ("absolute_return_tail_ratio", "absolute-tail ratio"),
    ("taylor_gap", "Taylor gap"),
    ("p_cap_exceedance_fraction", "P cap fraction"),
)


class FactorialMarkdownReport:
    """Render results without overstating illustrative acceptance bands."""

    def build(
        self,
        result: FactorialExperimentResult,
        latent_calibration: Mapping[str, object] | None = None,
    ) -> str:
        """Create a self-contained experiment ledger and comparison report."""
        manifest = result.manifest
        lines = [
            f"# {manifest.get('study_name', 'M2/M3/M2+O/M4 factorial pilot')}",
            "",
            "## Status and interpretation boundary",
            "",
            "This is a physical-measure capacity study on synthetic model paths. The "
            "acceptance bands are declared proxies, not empirical SPX targets. "
            "Accordingly, acceptance rates compare the four cells under one declared "
            "sampling measure; they do not yet validate the model against market data.",
            "",
            self._design_summary(manifest),
        ]
        if latent_calibration is not None:
            lines.extend(("", self._latent_summary(latent_calibration)))
        lines.extend(
            (
                "",
                "## Acceptance by independent cache seed",
                "",
                self._acceptance_table(result.seed_summary),
                "",
                "## Acceptance by structural scenario",
                "",
                self._scenario_table(result.scenario_summary),
                "",
                "## Stylised-fact comparison",
                "",
                "Entries are fixed-design means across seeds; the value after `±` is the "
                "between-seed standard deviation.",
                "",
                self._metric_table(result.metric_summary),
                "",
                "## Individual gate pass rates",
                "",
                "These rates identify the bottlenecks hidden by joint acceptance.",
                "",
                self._criterion_table(result.criterion_summary),
                "",
                "## Post-hoc threshold sensitivity",
                "",
                "These deliberately stricter thresholds are diagnostics, not preregistered "
                "acceptance criteria or empirical SPX estimates. They test whether the "
                "declared proxy gates are masking weak magnitudes.",
                "",
                self._stress_table(self.stress_summary(result.runs)),
                "",
                "## Exploratory M2 feedback-curvature screen",
                "",
                "This post-hoc quartile screen is descriptive, not causal. It is included "
                "to choose the next parameter box rather than to validate a fitted effect.",
                "",
                self._curvature_table(self.feedback_curvature_screen(result.runs)),
                "",
                "## Paired optional-component effects",
                "",
                "Contrasts use identical primitive paths, structural scenarios, and core "
                "readout draws. Positive acceptance effects are improvements; for distance, "
                "negative effects are improvements. Other signs retain the metric's meaning.",
                "",
                self._effect_table(result.effect_summary),
                "",
                "## Evidence summary",
                "",
                self._pilot_observations(result),
                "",
                "## Files and continuation rule",
                "",
                "Raw rows, per-seed checkpoints, seed/scenario/gate summaries, metric "
                "summaries, paired effects, and the complete manifest are stored beside this "
                "report. Keep an "
                "optional component only if its improvement exceeds between-seed uncertainty "
                "and does not materially degrade leverage, Zumbach, return ACF, or cap behaviour.",
                "",
                f"With only {len(manifest['seeds'])} seeds, uncertainty intervals are "
                "necessarily wide. A later market-data study must replace the proxy bands "
                "and retain an untouched validation period.",
            )
        )
        return "\n".join(lines).rstrip() + "\n"

    def save(
        self,
        path: str | Path,
        result: FactorialExperimentResult,
        latent_calibration: Mapping[str, object] | None = None,
    ) -> Path:
        """Write the report after all aggregate tables have been validated."""
        output = Path(path)
        output.write_text(self.build(result, latent_calibration), encoding="utf-8")
        return output

    @staticmethod
    def _design_summary(manifest: Mapping[str, object]) -> str:
        """State the fixed design and replication dimensions explicitly."""
        simulation = manifest["simulation"]
        return "\n".join(
            (
                "## Experimental design",
                "",
                f"- Cache seeds: `{manifest['seeds']}`.",
                f"- Full readout draws: `{manifest['candidate_count']}` using "
                f"`{manifest['sampler']}` with design seed `{manifest['design_seed']}`.",
                f"- Parameter space: `{manifest.get('parameter_space', {}).get('name', 'unrecorded')}`.",
                f"- Structural scenarios: `{len(manifest['structural_scenarios'])}`.",
                f"- Simulation: `{simulation}`.",
                "- Cells: M2 (neither optional loading), M3 (echo only), M2+O "
                "(orthogonal energy only), and M4 (both).",
                "- The echo consumes only the common feedback, spike, and clustering banks, "
                "so the orthogonal switch cannot silently redefine M3.",
            )
        )

    @staticmethod
    def _latent_summary(calibration: Mapping[str, object]) -> str:
        """Record the architecture-only fit separately from full-model evidence."""
        return "\n".join(
            (
                "## Latent architecture calibration",
                "",
                f"- Optimizer seeds: `{calibration.get('optimizer_seeds', [])}`; selected "
                f"seed: `{calibration.get('selected_seed')}`.",
                f"- Objective: `{_format_number(calibration.get('objective'))}`.",
                f"- Objective range: `{calibration.get('objective_range')}`.",
                f"- Optimizer status: `{calibration.get('optimizer_message')}`.",
                f"- Feedback rates: `{calibration.get('feedback_rate_range')}`.",
                f"- Clustering rates: `{calibration.get('cluster_rate_range')}`.",
                f"- Search-bound hits: `{calibration.get('boundary_hits', [])}`.",
                "",
                "This fit targets only the finite-band rough variogram and latent clustering "
                "curve. It does not determine leverage, Zumbach, tails, or Taylor effects. "
                "Endpoint hits require a wider-box or regularised follow-up before treating "
                "the architecture as identified.",
            )
        )

    @staticmethod
    def _acceptance_table(seed_summary: pd.DataFrame) -> str:
        """Show every seed rather than hiding replication variation."""
        pivot = seed_summary.pivot(index="seed", columns="variant", values="acceptance_rate")
        rows = []
        for seed, values in pivot.iterrows():
            row = [str(seed)]
            row.extend(_format_percent(values.get(item.value)) for item in VARIANT_ORDER)
            rows.append(row)
        headers = ["seed", *(item.value for item in VARIANT_ORDER)]
        return _markdown_table(headers, rows)

    @staticmethod
    def _metric_table(summary: pd.DataFrame) -> str:
        """Compact the principal diagnostics into one model-by-metric table."""
        rows = []
        for variant in (item.value for item in VARIANT_ORDER):
            row = [variant]
            subset = summary[summary.variant == variant].set_index("metric")
            for metric, _ in DISPLAY_METRICS:
                values = subset.loc[metric]
                row.append(_format_mean_sd(values["mean"], values["between_seed_sd"]))
            rows.append(row)
        headers = ["model", *(label for _, label in DISPLAY_METRICS)]
        return _markdown_table(headers, rows)

    @staticmethod
    def _scenario_table(summary: pd.DataFrame) -> str:
        """Expose structural sensitivity instead of averaging it away."""
        acceptance = summary[summary.metric == "acceptance_rate"]
        rows = []
        for scenario in sorted(acceptance.structural_scenario_id.unique()):
            row = [str(int(scenario))]
            subset = acceptance[acceptance.structural_scenario_id == scenario].set_index(
                "variant"
            )
            for variant in (item.value for item in VARIANT_ORDER):
                values = subset.loc[variant]
                row.append(
                    _format_percent_mean_sd(values["mean"], values["between_seed_sd"])
                )
            rows.append(row)
        return _markdown_table(["scenario", *(item.value for item in VARIANT_ORDER)], rows)

    @staticmethod
    def _criterion_table(summary: pd.DataFrame) -> str:
        """Report gate pass rates averaged only after seed-level evaluation."""
        criteria = tuple(summary.criterion.drop_duplicates())
        rows = []
        for variant in (item.value for item in VARIANT_ORDER):
            subset = summary[summary.variant == variant].set_index("criterion")
            row = [variant]
            for criterion in criteria:
                values = subset.loc[criterion]
                row.append(
                    _format_percent_mean_sd(values["mean"], values["between_seed_sd"])
                )
            rows.append(row)
        return _markdown_table(["model", *criteria], rows)

    @staticmethod
    def stress_summary(runs: pd.DataFrame) -> pd.DataFrame:
        """Evaluate stricter diagnostic thresholds separately within every seed."""
        seed_rows = []
        for (seed, variant), frame in runs.groupby(["seed", "variant"], sort=False):
            masks = FactorialMarkdownReport._stress_masks(frame)
            seed_rows.extend(
                {
                    "seed": int(seed),
                    "variant": variant,
                    "check": name,
                    "pass_rate": float(mask.mean()),
                }
                for name, mask in masks.items()
            )
        per_seed = pd.DataFrame(seed_rows)
        rows = []
        for (variant, check), frame in per_seed.groupby(["variant", "check"], sort=False):
            values = frame.pass_rate.to_numpy(dtype=float)
            rows.append(
                {
                    "variant": variant,
                    "check": check,
                    "mean": float(np.mean(values)),
                    "between_seed_sd": float(np.std(values, ddof=1)),
                    "seed_min": float(np.min(values)),
                    "seed_max": float(np.max(values)),
                }
            )
        return pd.DataFrame(rows)

    @staticmethod
    def _stress_masks(frame: pd.DataFrame) -> dict[str, pd.Series]:
        """Declare transparent magnitude checks for the first capacity pilot."""
        masks = {
            "rough H in [0.08, 0.12]": frame.rough_hurst.between(0.08, 0.12),
            "rank leverage < -0.01": frame.leverage_rank < -0.01,
            "rank Zumbach > 0.02": frame.zumbach_rank > 0.02,
            "rank Zumbach > 0.05": frame.zumbach_rank > 0.05,
            "excess kurtosis > 3": frame.excess_kurtosis > 3.0,
            "excess kurtosis > 5": frame.excess_kurtosis > 5.0,
        }
        masks["joint H/leverage/Z>0.02/kurtosis>3"] = (
            masks["rough H in [0.08, 0.12]"]
            & masks["rank leverage < -0.01"]
            & masks["rank Zumbach > 0.02"]
            & masks["excess kurtosis > 3"]
        )
        return masks

    @staticmethod
    def _stress_table(summary: pd.DataFrame) -> str:
        """Display seed-averaged stress-check rates without hiding dispersion."""
        checks = tuple(summary.check.drop_duplicates())
        rows = []
        for variant in (item.value for item in VARIANT_ORDER):
            subset = summary[summary.variant == variant].set_index("check")
            row = [variant]
            for check in checks:
                values = subset.loc[check]
                row.append(
                    _format_percent_mean_sd(values["mean"], values["between_seed_sd"])
                )
            rows.append(row)
        return _markdown_table(["model", *checks], rows)

    @staticmethod
    def feedback_curvature_screen(runs: pd.DataFrame) -> pd.DataFrame:
        """Summarise M2 by fixed-design feedback-curvature quantiles and cache seed."""
        frame = runs[runs.variant == "M2"].copy()
        quantile_count = min(4, int(frame.feedback_curvature.nunique()))
        frame["curvature_bin"] = pd.qcut(
            frame.feedback_curvature,
            q=quantile_count,
            duplicates="drop",
        ).astype(str)
        metrics = (
            "rough_hurst",
            "excess_kurtosis",
            "leverage_rank",
            "zumbach_rank",
            "distance",
            "accepted",
        )
        per_seed = frame.groupby(
            ["seed", "curvature_bin"], as_index=False, observed=True
        )[list(metrics)].mean()
        ranges = frame.groupby("curvature_bin", observed=True).feedback_curvature.agg(
            ["min", "max"]
        )
        rows = []
        for curvature_bin, group in per_seed.groupby("curvature_bin", sort=False):
            row: dict[str, object] = {
                "curvature_bin": curvature_bin,
                "feedback_curvature_min": float(ranges.loc[curvature_bin, "min"]),
                "feedback_curvature_max": float(ranges.loc[curvature_bin, "max"]),
            }
            for metric in metrics:
                values = group[metric].to_numpy(dtype=float)
                row[f"{metric}_mean"] = float(np.mean(values))
                row[f"{metric}_between_seed_sd"] = float(np.std(values, ddof=1))
            rows.append(row)
        return pd.DataFrame(rows).sort_values("feedback_curvature_min").reset_index(drop=True)

    @staticmethod
    def _curvature_table(summary: pd.DataFrame) -> str:
        """Display the principal M2 trade-offs across curvature quartiles."""
        metrics = (
            "rough_hurst",
            "excess_kurtosis",
            "leverage_rank",
            "zumbach_rank",
            "accepted",
        )
        rows = []
        for _, values in summary.iterrows():
            row = [
                f"[{values['feedback_curvature_min']:.4g}, "
                f"{values['feedback_curvature_max']:.4g}]"
            ]
            row.extend(
                _format_mean_sd(
                    values[f"{metric}_mean"],
                    values[f"{metric}_between_seed_sd"],
                )
                for metric in metrics
            )
            rows.append(row)
        return _markdown_table(["feedback curvature", *metrics], rows)

    @staticmethod
    def _effect_table(summary: pd.DataFrame) -> str:
        """Report selected raw-scale factorial contrasts with seed uncertainty."""
        selected = (
            "acceptance_rate",
            "distance",
            "rough_hurst",
            "leverage_rank",
            "zumbach_rank",
            "excess_kurtosis",
        )
        rows = []
        for effect in summary.effect.drop_duplicates():
            subset = summary[summary.effect == effect].set_index("metric")
            row = [effect]
            for metric in selected:
                values = subset.loc[metric]
                row.append(_format_mean_sd(values["mean"], values["between_seed_sd"]))
            rows.append(row)
        return _markdown_table(["effect", *selected], rows)

    @staticmethod
    def _pilot_observations(result: FactorialExperimentResult) -> str:
        """Turn principal paired estimates into cautious, auditable statements."""
        metrics = result.metric_summary.set_index(["variant", "metric"])
        effects = result.effect_summary.set_index(["effect", "metric"])
        criteria = result.criterion_summary.set_index(["variant", "criterion"])
        acceptance = {
            variant: metrics.loc[(variant, "acceptance_rate"), "mean"]
            for variant in (item.value for item in VARIANT_ORDER)
        }
        echo = effects.loc[("echo_at_O0", "distance")]
        echo_leverage = effects.loc[("echo_at_O0", "leverage_rank")]
        echo_zumbach = effects.loc[("echo_at_O0", "zumbach_rank")]
        orthogonal = effects.loc[("orthogonal_at_echo0", "distance")]
        zumbach_gate = criteria.loc[("M2", "zumbach_rank")]
        zumbach_spread = zumbach_gate["seed_max"] - zumbach_gate["seed_min"]
        zumbach_interpretation = (
            "This is a warning against relying on one simulation seed."
            if zumbach_spread > 0.20
            else "The sign gate is comparatively stable across the declared seeds."
        )
        accepted = result.runs.pivot(
            index=["seed", "structural_scenario_id", "candidate_id"],
            columns="variant",
            values="accepted",
        ).astype(bool)
        echo_gained = int(((~accepted["M2"]) & accepted["M3"]).sum())
        echo_lost = int((accepted["M2"] & (~accepted["M3"])).sum())
        orthogonal_gained = int(((~accepted["M2"]) & accepted["M2+O"]).sum())
        orthogonal_lost = int((accepted["M2"] & (~accepted["M2+O"])).sum())
        stress = FactorialMarkdownReport.stress_summary(result.runs).set_index(
            ["variant", "check"]
        )
        tight_roughness = stress.loc[("M2", "rough H in [0.08, 0.12]"), "mean"]
        strong_tails = stress.loc[("M2", "excess kurtosis > 5"), "mean"]
        strict_joint = stress.loc[
            ("M2", "joint H/leverage/Z>0.02/kurtosis>3"), "mean"
        ]
        curvature = FactorialMarkdownReport.feedback_curvature_screen(result.runs)
        best_curvature = curvature.loc[curvature.distance_mean.idxmin()]
        return "\n".join(
            (
                f"- Joint acceptance was M2 `{_format_percent(acceptance['M2'])}`, M3 "
                f"`{_format_percent(acceptance['M3'])}`, M2+O "
                f"`{_format_percent(acceptance['M2+O'])}`, and M4 "
                f"`{_format_percent(acceptance['M4'])}` under the declared proxy bands.",
                f"- The echo-only contrast changed distance by `{_format_ci(echo)}`. It "
                f"changed rank leverage by `{_format_ci(echo_leverage)}` and rank Zumbach "
                f"by `{_format_ci(echo_zumbach)}`; signs must be judged against the desired "
                "negative leverage and positive Zumbach directions.",
                f"- The orthogonal-only distance contrast was `{_format_ci(orthogonal)}`. "
                "An effect of this size should not justify the additional bank unless it "
                "becomes stable under a broader loading range or a targeted conditional test.",
                f"- At the row level, the echo changed `{echo_gained}` failures into passes "
                f"and `{echo_lost}` passes into failures. The orthogonal block changed "
                f"`{orthogonal_gained}` failures into passes and `{orthogonal_lost}` passes "
                "into failures. This distinguishes a small net rate change from uniformly "
                "better candidates.",
                f"- For M2, the Zumbach gate pass rate averaged "
                f"`{_format_percent(zumbach_gate['mean'])}` and ranged from "
                f"`{_format_percent(zumbach_gate['seed_min'])}` to "
                f"`{_format_percent(zumbach_gate['seed_max'])}` across seeds. "
                f"{zumbach_interpretation}",
                f"- The tighter post-hoc checks are substantially less favourable: M2 "
                f"reached H in [0.08, 0.12] in `{_format_percent(tight_roughness)}` of "
                f"rows, excess kurtosis above 5 in `{_format_percent(strong_tails)}`, and "
                f"the stricter joint check in `{_format_percent(strict_joint)}`. Therefore "
                "the proxy joint rate must not be read as evidence that all "
                "target magnitudes have already been obtained.",
                f"- The lowest-distance M2 feedback-curvature bin is the clearest next search "
                f"region: mean H was `{_format_number(best_curvature['rough_hurst_mean'])}`, "
                f"excess kurtosis `{_format_number(best_curvature['excess_kurtosis_mean'])}`, "
                f"rank leverage `{_format_number(best_curvature['leverage_rank_mean'])}`, "
                f"and joint proxy acceptance "
                f"`{_format_percent(best_curvature['accepted_mean'])}`. Rank Zumbach was "
                f"`{_format_number(best_curvature['zumbach_rank_mean'])}`, so this is a "
                "promising narrowing direction, not a complete solution.",
            )
        )


def _markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    """Render a simple CommonMark table without an additional dependency."""
    escaped_headers = [_escape_cell(item) for item in headers]
    lines = [
        "| " + " | ".join(escaped_headers) + " |",
        "| " + " | ".join("---" for _ in escaped_headers) + " |",
    ]
    lines.extend("| " + " | ".join(_escape_cell(item) for item in row) + " |" for row in rows)
    return "\n".join(lines)


def _escape_cell(value: object) -> str:
    """Protect Markdown separators inside compact table cells."""
    return str(value).replace("|", "\\|").replace("\n", " ")


def _format_number(value: object) -> str:
    """Use compact scientific formatting for heterogeneous scalar results."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return f"{number:.4g}" if np.isfinite(number) else "NA"


def _format_mean_sd(mean: object, standard_deviation: object) -> str:
    """Display the estimate and visible between-seed variation."""
    return f"{_format_number(mean)} ± {_format_number(standard_deviation)}"


def _format_percent(value: object) -> str:
    """Format an acceptance fraction as a percentage."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    return f"{100.0 * number:.1f}%" if np.isfinite(number) else "NA"


def _format_percent_mean_sd(mean: object, standard_deviation: object) -> str:
    """Display a rate and its between-seed dispersion in percentage points."""
    try:
        center, spread = float(mean), float(standard_deviation)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(center) or not np.isfinite(spread):
        return "NA"
    return f"{100.0 * center:.1f}% ± {100.0 * spread:.1f}pp"


def _format_ci(values: pd.Series) -> str:
    """Display a seed-mean effect and its small-sample confidence interval."""
    mean = _format_number(values["mean"])
    lower = _format_number(values["ci_lower"])
    upper = _format_number(values["ci_upper"])
    return f"{mean} (95% seed CI [{lower}, {upper}])"
