"""Make observation/termination changes to a pipeline-generated GH input.

No initial-data, gauge, damping, domain, spatial-discretization or control-system
parameters are changed. The resulting input MUST pass the actual C++ parser.
"""
import argparse
import json
from pathlib import Path

import yaml


def prepare(source, destination):
    metadata, data = list(yaml.safe_load_all(source.read_text()))
    changes = []

    def replace(parent, key, value, path):
        changes.append({"path": path, "old": parent.get(key), "new": value})
        parent[key] = value

    replace(data["Evolution"], "InitialSlabSize", 1e-4, "Evolution.InitialSlabSize")
    replace(data["Evolution"], "InitialTimeStep", 1e-5, "Evolution.InitialTimeStep")
    replace(data["Observers"], "VolumeFileName", "GhThreeTimeVolume", "Observers.VolumeFileName")
    replace(data["Observers"], "ReductionFileName", "GhThreeTimeReductions", "Observers.ReductionFileName")
    replace(data["Observers"], "SurfaceFileName", "GhThreeTimeSurfaces", "Observers.SurfaceFileName")
    fields = ["SpacetimeMetric", "Pi", "Phi", "GaugeH", "SpacetimeDerivGaugeH",
              "Lapse", "Shift", "SpatialMetric", "GaugeConstraint"]
    events = [
        {"Trigger": "Always", "Events": [{"ObserveFields": {
            "SubfileName": "Validation", "VariablesToObserve": fields,
            "InterpolateToMesh": "None", "ProjectToMesh": "None",
            "CoordinatesFloatingPointType": "Double", "FloatingPointTypes": ["Double"],
            "BlocksToObserve": "All"}}]},
        {"Trigger": "Always", "Events": [{"ObserveNorms": {
            "SubfileName": "Constraints", "TensorsToObserve": [
                {"Name": "Lapse", "NormType": "Min", "Components": "Individual"},
                {"Name": "PointwiseL2Norm(GaugeConstraint)", "NormType": "L2Norm", "Components": "Sum"},
                {"Name": "PointwiseL2Norm(ThreeIndexConstraint)", "NormType": "L2Norm", "Components": "Sum"}],
            "BlocksToObserve": "All", "ReductionScope": "WholeDomain"}}]},
        {"Trigger": {"TimeCompares": {"Comparison": "GreaterThanOrEqualTo", "Value": 2e-4}},
         "Events": ["Completion"]},
    ]
    # This also removes CleanFunctionsOfTime to retain the full short history.
    replace(data, "EventsAndTriggersAtSlabs", events, "EventsAndTriggersAtSlabs")
    # Disable only the user-configured CCE extraction event, not control systems.
    replace(data, "EventsAndDenseTriggers", [], "EventsAndDenseTriggers")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(yaml.safe_dump_all([metadata, data], sort_keys=False))
    (destination.parent / "input-changes.json").write_text(json.dumps(changes, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    prepare(args.source, args.destination)
