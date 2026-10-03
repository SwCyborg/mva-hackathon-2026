#!/usr/bin/env python3
"""
Validador local del submission Track 1 — replica la logica de evaluation.py oficial.
Verifica que el CSV parsea sin errores y simula el scoring contra los escenarios probables.
NO es el evaluador oficial (no tengo el answer key), es un chequeo de integridad.
"""
import csv, sys

RANK_POINT_TIERS = [(1, 100), (3, 50), (5, 25), (10, 10)]

def _parse_variant(chrom, pos, ref, alt):
    if not chrom or not pos:
        return None
    return (str(chrom).strip(), int(pos), str(ref).strip().upper(), str(alt).strip().upper())

def load_submission(path):
    """Replica load_submission del evaluador. Valida reglas."""
    by_proband = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            pid = r["proband_id"].strip()
            v1 = _parse_variant(r["chrom_1"], r["pos_1"], r["ref_1"], r["alt_1"])
            if v1 is None:
                raise ValueError(f"Fila sin variante primaria: {r}")
            variants = {v1}
            v2 = _parse_variant(r["chrom_2"], r["pos_2"], r["ref_2"], r["alt_2"])
            if v2 is not None:
                variants.add(v2)
            epcr = float(r["epcr"])
            if not (0 < epcr <= 1):
                raise ValueError(f"epcr fuera de (0,1]: {epcr}")
            ft = r.get("finding_type", "primary").strip()
            if ft not in ("primary", "secondary"):
                raise ValueError(f"finding_type invalido: {ft}")
            by_proband.setdefault(pid, []).append((frozenset(variants), epcr))
    for pid, rows in by_proband.items():
        if len(rows) > 10:
            raise ValueError(f"{pid}: {len(rows)} filas > 10")
    # ordenar desc por epcr
    for pid in by_proband:
        by_proband[pid].sort(key=lambda r: -r[1])
    return by_proband

def _rank_to_points(rank):
    for t, p in RANK_POINT_TIERS:
        if rank <= t:
            return p
    return 0

def score_proband(true_variants, rows):
    ts = frozenset(true_variants)
    is_ch = len(ts) == 2
    fr = pr = None
    for i, (v, e) in enumerate(rows, 1):
        if v == ts:
            fr = i; break
    if fr is None and is_ch:
        for i, (v, e) in enumerate(rows, 1):
            if len(v & ts) == 1:
                pr = i; break
    rp = _rank_to_points(fr) if fr else (0.5 * _rank_to_points(pr) if pr else 0.0)
    best = 0.0
    for t in sorted({r[1] for r in rows}, reverse=True):
        pred = set()
        for v, e in rows:
            if e >= t:
                pred |= set(v)
        tp = len(pred & ts); fp = len(pred - ts); fn = len(ts - pred)
        p = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        best = max(best, 2 * p * rec / (p + rec) if p + rec else 0.0)
    return rp, best, fr, pr

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "/home/sw-cyborg/mva_hackathon/Swcyborg_track1_submission.csv"
    print(f"=== VALIDANDO {path} ===")
    try:
        subs = load_submission(path)
    except Exception as e:
        print(f"🔴 FALLO DE VALIDACION: {e}")
        sys.exit(1)
    print(f"✅ Parsea OK. Probandos: {list(subs.keys())}")
    for pid, rows in subs.items():
        print(f"   {pid}: {len(rows)} filas")
        for i, (v, e) in enumerate(rows, 1):
            print(f"      rank{i} epcr={e:.2f} variants={sorted(v)}")

    # Simular escenarios de answer key
    V1=("chr15",40209701,"T","G"); V2=("chr15",40220612,"T","G")
    V1n=("15",40209701,"T","G");   V2n=("15",40220612,"T","G")
    print("\n=== SIMULACION DE SCORING (escenarios) ===")
    for name, true in [("key=par chr15",[V1,V2]),("key=par '15'",[V1n,V2n]),
                       ("key=V1 chr15",[V1]),("key=V1 '15'",[V1n]),("key=V2 chr15",[V2])]:
        rows = subs["PROBAND01"]
        rp, f, fr, pr = score_proband(true, rows)
        print(f"   {name:14s} -> rank_pts={rp:5.1f}  f_max={f:.3f}  full@{fr}  partial@{pr}")
    print("\n✅ VALIDACION COMPLETA")
