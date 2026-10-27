"""Print a compact table from an evaluate.py --compare output file (concatenated JSON objects)."""
import json, sys

txt = open(sys.argv[1]).read()
dec, i, objs = json.JSONDecoder(), 0, []
while (j := txt.find("{", i)) >= 0:
    try:
        o, e = dec.raw_decode(txt[j:]); objs.append(o); i = j + e
    except json.JSONDecodeError:
        i = j + 1
mets = {o["run_id"]: o for o in objs if "run_id" in o and "OR@0.1%" in o}
print(f"{'candidate':28s} {'OR@0.1%':>8s} {'dlogOR':>7s} {'fold1':>6s} {'fold2':>6s} {'jk_lb90':>7s} {'bt_lb90':>7s} {'dCMH':>6s} adopt")
for o in objs:
    if "candidate" not in o:
        continue
    m = mets.get(o["candidate"], {})
    f = o["delta_log_or_folds"]
    print(f"{o['candidate']:28s} {m.get('OR@0.1%', float('nan')):8.2f} {o['delta_log_or']:7.3f} {f['fold1']:6.2f} {f['fold2']:6.2f} "
          f"{o['paired_jackknife']['lb90_one_sided']:7.3f} {o['boot_delta_log_or_lb90']:7.3f} {o['delta_log_cmh']:6.3f} {o['adopt_rule_8_2']}")
