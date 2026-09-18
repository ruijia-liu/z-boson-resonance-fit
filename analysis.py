"""Reproducible educational Z-resonance fit; see docs/METHODS.md."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.constants import hbar, electron_volt
from scipy.optimize import curve_fit

HBAR_GEV_S = hbar / (1e9 * electron_volt)
ROOT = Path(__file__).resolve().parent


def cross_section(energy, mass, width, pole_cross_section):
    """Constant-width shape, normalized to sigma(E=mass), in nb.

    Energies, mass parameter and total width are in GeV (natural units).
    The E**2 numerator preserves the original coursework line-shape choice.
    """
    energy = np.asarray(energy, dtype=float)
    return pole_cross_section * energy**2 * width**2 / (
        (energy**2 - mass**2)**2 + mass**2 * width**2)


def lifetime(width):
    if not np.isfinite(width) or width <= 0:
        raise ValueError("Width must be finite and positive")
    return HBAR_GEV_S / width


def load_data(paths):
    """Retain source row numbers and audit every exclusion; never deduplicate observations."""
    accepted, audit = [], []
    for path in paths:
        with Path(path).open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(stream)
            header = next(reader, None)
            if header is None or len(header) != 3:
                raise ValueError(f"Expected a three-column CSV header: {Path(path).name}")
            for line, row in enumerate(reader, 2):
                reason = ""
                try:
                    if len(row) != 3:
                        raise ValueError("Expected three fields")
                    values = np.array([float(v) for v in row])
                    if not np.isfinite(values).all():
                        reason = "nonfinite"
                    elif values[0] <= 0 or values[2] <= 0:
                        reason = "nonpositive_energy_or_uncertainty"
                    elif not 85 < values[0] < 95:
                        reason = "outside_coursework_energy_window"
                    elif not 0 < values[1] < 3:
                        reason = "outside_coursework_cross_section_window"
                except ValueError:
                    reason = "nonnumeric_or_wrong_column_count"
                record = {"source": Path(path).name, "line": line,
                          "raw": row, "status": reason or "candidate"}
                audit.append(record)
                if not reason:
                    accepted.append(values)
    data = np.array(accepted, dtype=float).reshape(-1, 3)
    return data, audit


def fit(data):
    data = np.asarray(data, dtype=float)
    if (data.ndim != 2 or data.shape[1] != 3 or len(data) <= 3
            or not np.isfinite(data).all() or np.any(data[:, 2] <= 0)
            or len(np.unique(data[:, 0])) < 4):
        raise ValueError("Need more than three valid observations at distinct energies")
    parameters, covariance = curve_fit(
        cross_section, data[:, 0], data[:, 1], sigma=data[:, 2],
        p0=[91.2, 2.5, 2.0], absolute_sigma=True,
        bounds=([85.0, 0.01, 0.0], [95.0, 20.0, np.inf]), maxfev=20000)
    if not np.isfinite(covariance).all():
        raise ValueError("Fit covariance is not finite")
    residuals = (data[:, 1] - cross_section(data[:, 0], *parameters)) / data[:, 2]
    errors = np.sqrt(np.diag(covariance))
    tau = lifetime(parameters[1])
    return {"parameters": dict(zip(["mass_GeV", "width_GeV", "pole_cross_section_nb"], parameters.tolist())),
            "standard_errors": dict(zip(["mass_GeV", "width_GeV", "pole_cross_section_nb"], errors.tolist())),
            "covariance": covariance.tolist(), "chi_squared": float(residuals @ residuals),
            "degrees_of_freedom": len(data) - 3,
            "reduced_chi_squared": float(residuals @ residuals / (len(data) - 3)),
            "lifetime_s": tau, "lifetime_standard_error_s": tau * errors[1] / parameters[1],
            "observations": len(data)}


def parameters(result):
    return np.array(list(result["parameters"].values()))


def analyze(data, audit, clip_sigma=5.0):
    initial = fit(data)
    residuals = (data[:, 1] - cross_section(data[:, 0], *parameters(initial))) / data[:, 2]
    keep = np.ones(len(data), dtype=bool) if clip_sigma is None else np.abs(residuals) < clip_sigma
    for record, retain, residual in zip((r for r in audit if r["status"] == "candidate"), keep, residuals):
        record["status"] = "retained" if retain else "prefit_residual_cut"
        record["prefit_standardized_residual"] = float(residual)
    final = fit(data[keep])
    return initial, final, keep


def profile_grid(data, masses, widths):
    """Analytically optimize the nonnegative linear amplitude at every grid point."""
    mass, width = np.meshgrid(masses, widths)
    shape = cross_section(data[:, 0, None, None], mass, width, 1.0)
    weights = 1 / data[:, 2, None, None]**2
    amplitude = np.maximum(0, (weights * shape * data[:, 1, None, None]).sum(0)
                           / (weights * shape**2).sum(0))
    chi2 = (((data[:, 1, None, None] - amplitude * shape) / data[:, 2, None, None])**2).sum(0)
    return mass, width, chi2


def plot_results(data, keep, result, output):
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    p = parameters(result)
    fig, (ax, residual_ax) = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                                         gridspec_kw={"height_ratios": [3, 1]}, layout="constrained")
    ax.errorbar(*data[keep, :2].T, yerr=data[keep, 2], fmt="o", ms=3.5,
                color="#236c9c", alpha=.85, label="Retained coursework observations")
    if (~keep).any():
        ax.errorbar(*data[~keep, :2].T, yerr=data[~keep, 2], fmt="x", color="#c55732",
                    label="Excluded by one-pass residual cut")
    energy = np.linspace(data[:, 0].min(), data[:, 0].max(), 700)
    ax.plot(energy, cross_section(energy, *p), color="#172f45", lw=2, label="Fitted resonance")
    ax.set(ylabel="Cross section (nb)", title="Z resonance | coursework data")
    ax.legend(fontsize=9)
    residuals = (data[:, 1] - cross_section(data[:, 0], *p)) / data[:, 2]
    residual_ax.scatter(data[keep, 0], residuals[keep], s=13, color="#236c9c")
    residual_ax.scatter(data[~keep, 0], residuals[~keep], s=25, marker="x", color="#c55732")
    residual_ax.axhline(0, color="black", lw=.8)
    residual_ax.axhline(2, color="gray", ls="--", lw=.8)
    residual_ax.axhline(-2, color="gray", ls="--", lw=.8)
    residual_ax.set(xlabel="Centre-of-mass energy (GeV)", ylabel="Residual / error")
    fig.savefig(output / "resonance_fit.png", dpi=180)
    plt.close(fig)
    errors = np.array(list(result["standard_errors"].values()))
    masses = np.linspace(p[0] - 4*errors[0], p[0] + 4*errors[0], 161)
    widths = np.linspace(max(.001, p[1] - 4*errors[1]), p[1] + 4*errors[1], 161)
    mass, width, chi2 = profile_grid(data[keep], masses, widths)
    fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
    contours = ax.contour(mass, width, chi2 - result["chi_squared"], levels=[2.30, 6.18],
                          colors=["#236c9c", "#c55732"])
    ax.clabel(contours, fmt={2.30: "68.3% nominal", 6.18: "95.4% nominal"}, fontsize=9)
    ax.scatter(p[0], p[1], color="#172f45", s=30, label="Best fit")
    ax.set(xlabel="Mass parameter (GeV)", ylabel="Total width (GeV)",
           title="Profile chi-square | amplitude re-fitted")
    ax.legend()
    fig.savefig(output / "profile_chi_squared.png", dpi=180)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", nargs="+", type=Path,
                        default=[ROOT / "data" / f"z_boson_data_{n}.csv" for n in (1, 2)])
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    parser.add_argument("--no-clip", action="store_true", help="Disable residual clipping; input windows still apply")
    args = parser.parse_args(argv)
    data, audit = load_data(args.data)
    initial, final, keep = analyze(data, audit, None if args.no_clip else 5.0)
    args.output.mkdir(parents=True, exist_ok=True)
    plot_results(data, keep, final, args.output)
    report = {"dataset": "Supplied PHYS10362 coursework CSVs; experimental provenance unverified",
              "model": "Constant-width, E-squared-numerator shape; free pole cross section",
              "residual_cut_sigma": None if args.no_clip else 5,
              "before_residual_cut": initial, "final": final,
              "input_rows": len(audit), "input_window_and_validity_exclusions": len(audit)-len(data),
              "residual_cut_exclusions": int((~keep).sum()),
              "inputs": [{"name": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in args.data],
              "versions": {"python": platform.python_version(), "numpy": np.__version__,
                           "scipy": scipy.__version__, "matplotlib": matplotlib.__version__}}
    (args.output / "fit_results.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    with (args.output / "row_audit.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["source", "line", "energy_raw", "cross_section_raw", "uncertainty_raw", "status", "prefit_standardized_residual"])
        for row in audit:
            writer.writerow([row["source"], row["line"], *(row["raw"]+[""]*3)[:3], row["status"], row.get("prefit_standardized_residual", "")])
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
