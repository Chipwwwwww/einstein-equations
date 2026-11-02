"""Inspect this XCTS output smoke test. All diagnostics are NUMERICAL.

Observation values are elliptic AMR iteration labels, NEVER evolution times.
Optional map checks require bindings from the pinned SpECTRE distribution.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def matrix(group, prefix):
    n = len(group["Lapse"])
    result = np.empty((n, 3, 3))
    for i, a in enumerate("xyz"):
        for j, b in enumerate("xyz"):
            suffix = a + b if i >= j else b + a
            result[:, i, j] = group[prefix + "_" + suffix][:]
    return result


def inspect(run_dir, check_maps=False):
    run_dir = Path(run_dir)
    result = {"classification": "EXECUTED / NUMERICAL",
              "physical_evolution_interval": None,
              "common_horizon_formation_tested": False,
              "files": {}, "observations": []}
    for path in sorted(run_dir.glob("*.h5")):
        result["files"][path.name] = {
            "bytes": path.stat().st_size, "sha256": sha256(path)}
    path = run_dir / "BbhVolume0.h5"
    with h5py.File(path, "r") as file:
        vol = file["VolumeData.vol"]
        result["serialized_geometry"] = {
            k: {"shape": list(vol[k].shape), "dtype": str(vol[k].dtype)}
            for k in ("domain", "global_functions_of_time")}
        groups = sorted((g for g in vol.values() if isinstance(g, h5py.Group)),
                        key=lambda g: g.attrs["observation_value"])
        for group in groups:
            extents = group["total_extents"][:].reshape(-1, 3)
            lapse = group["Lapse"][:]
            gamma = matrix(group, "SpatialMetric")
            beta = np.column_stack([group["ShiftExcess_" + a][:] for a in "xyz"])
            metric = np.empty((len(lapse), 4, 4))
            metric[:, 1:, 1:] = gamma
            metric[:, 0, 1:] = metric[:, 1:, 0] = np.einsum("nij,nj->ni", gamma, beta)
            metric[:, 0, 0] = -lapse**2 + np.einsum("ni,ni->n", beta, metric[:, 0, 1:])
            eigen = np.linalg.eigvalsh(metric)
            result["observations"].append({
                "path": group.name,
                "AMR_iteration_label_NOT_time": float(group.attrs["observation_value"]),
                "elements": len(extents), "grid_points": int(np.prod(extents, axis=1).sum()),
                "extents_present": np.unique(extents, axis=0).tolist(),
                "datasets": {k: {"shape": list(v.shape), "dtype": str(v.dtype)}
                             for k, v in group.items() if isinstance(v, h5py.Dataset)},
                "fields_present": {p: any(k == p or k.startswith(p + "_") for k in group)
                                   for p in ("SpacetimeMetric", "Pi", "Phi", "GaugeH",
                                             "SpacetimeDerivGaugeH", "Lapse", "Shift",
                                             "SpatialMetric", "ExtrinsicCurvature")},
                "min_lapse": float(lapse.min()),
                "min_spatial_metric_eigenvalue": float(np.linalg.eigvalsh(gamma).min()),
                "reconstructed_metric_shift": "ShiftExcess (pipeline inertial-shift convention)",
                "nodes_with_signature_minus_plus_plus_plus": int(np.count_nonzero(
                    (eigen[:, 0] < 0) & np.all(eigen[:, 1:] > 0, axis=1))),
                "all_metric_values_finite": bool(np.isfinite(metric).all()),
                "max_abs_hamiltonian_constraint": float(np.max(np.abs(group["HamiltonianConstraint"][:]))),
            })
        header = str(file.attrs["header.hdr"][0])
        result["compiled_library_versions"] = header.split("############### library versions ###############", 1)[-1].strip()
    with h5py.File(run_dir / "BbhReductions.h5", "r") as file:
        result["reductions"] = {
            name: {"legend": list(file[name].attrs["Legend"]), "data": file[name][:].tolist()}
            for name in ("ConstraintDiagnostics.dat", "NewtonRaphsonResiduals.dat")}
    if (run_dir / "Horizons.h5").exists():
        with h5py.File(run_dir / "Horizons.h5", "r") as file:
            result["individual_horizons"] = {
                name: {"legend": list(file[name].attrs["Legend"]), "data": file[name][:].tolist()}
                for name in ("AhA.dat", "AhB.dat")}
    if check_maps:
        result["map_checks"] = inspect_maps(path)
    return result


def inspect_maps(path):
    import spectre.IO.H5 as sh5
    from spectre.IO.H5.IterElements import iter_elements
    from spectre.Domain import deserialize_domain, deserialize_functions_of_time

    def bounds(fot):
        return [float(x) if np.isfinite(x) else str(x) for x in fot.time_bounds()]

    with sh5.H5File(str(path), "r") as file:
        vol = file.get_vol("/VolumeData")
        domain = deserialize_domain[3](vol.get_domain())
        global_fots = deserialize_functions_of_time(vol.get_global_functions_of_time())
        result = {"blocks": len(domain.blocks), "block_names": [b.name for b in domain.blocks],
                  "global_function_bounds": {k: bounds(f) for k, f in global_fots.items()},
                  "observations": []}
        for obs in vol.list_observation_ids():
            value = vol.get_observation_value(obs)
            fots = deserialize_functions_of_time(vol.get_functions_of_time(obs))
            mismatch = velocity = 0.0
            for el, data in iter_elements(vol, obs, tensor_components=["InertialCoordinates_" + a for a in "xyz"]):
                mismatch = max(mismatch, float(np.max(np.abs(np.asarray(el.inertial_coordinates) - data))))
                x1 = np.asarray(el.map(el.logical_coordinates, .25, global_fots))
                x2 = np.asarray(el.map(el.logical_coordinates, .75, global_fots))
                velocity = max(velocity, float(np.max(np.abs((x2 - x1) / .5))))
            result["observations"].append({
                "AMR_iteration_label_NOT_time": value,
                "bases": sorted({str(b) for row in vol.get_bases(obs) for b in row}),
                "quadratures": sorted({str(q) for row in vol.get_quadratures(obs) for q in row}),
                "observation_function_bounds": {k: bounds(f) for k, f in fots.items()},
                "max_coordinate_reconstruction_difference": mismatch,
                "static_map_finite_difference_velocity_max_abs": velocity,
                "max_function_first_derivative_at_observation": max(
                    float(np.max(np.abs(np.asarray(f.func_and_deriv(value)[1])))) for f in fots.values()),
            })
        return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--spectre-maps", action="store_true")
    args = parser.parse_args()
    print(json.dumps(inspect(args.run_dir, args.spectre_maps), indent=2, allow_nan=False))
