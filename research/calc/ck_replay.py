"""Run: python -m research.calc.ck_replay waist ext. Replay the ext and waist certificates in a fresh process at the committed source and compare every
certified field with the stored records."""
import json, subprocess, sys, time
from research.calc import ck_bound as ck
from research.calc import gap3_interval as gi
commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
dirty = subprocess.run(["git", "status", "--porcelain", "research/calc"], capture_output=True, text=True).stdout.strip()
out = {"commit": commit, "research_calc_clean": dirty == "", "runs": {}}
boxes = ck.profile_boxes()
for mode in sys.argv[1:]:
    cfg = ck.MODES[mode]
    stored = json.load(open(f"research/calc/data/{cfg['out']}"))
    region = boxes if cfg["region"] == "tube" else [b for b in boxes if b[3] <= ck.WAIST_ZMAX]
    t = time.time()
    r = ck.certify(cfg["s0"], cfg["radii"], ck.D_LO, gi.D_HOLE, cell=ck.CELL, boxes=region, sharp=cfg["sharp"],
                   im_split=cfg["im_split"], workers=4)
    keys = ["ok", "cells", "C", "K0", "K0w", "T_bg", "B_bg", "T_min", "sup_u_weighted_bound", "c", "c_loss", "c_zero",
            "rows_C", "rows_L", "schedule", "psi_abs_lower", "gam_op_upper"]
    diffs = [k for k in keys if r[k] != stored[k]]
    out["runs"][mode] = {"identical_fields": [k for k in keys if k not in diffs], "differing_fields": diffs,
                         "profile_identical": [list(b) for b in region] == stored["profile"],
                         "T_min": r["T_min"], "seconds": time.time() - t}
    print(mode, out["runs"][mode], flush=True)
json.dump(out, open("research/calc/data/ck_bound_replay.json", "w"), indent=1)
