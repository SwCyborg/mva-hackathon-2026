#!/usr/bin/env python3
"""
Llamado de CNV exónico de BUB1B por ratio de cobertura.
Método: mediana por exón / mediana global del gen -> ratio ~1.0 diploide.
         deleción het -> ratio ~0.5 en bloque de exones contiguos.
Usa el BAM SANO: b_larga/bam/region_38_43M.bam (GRCh38, 4 lanes, 54.7x parejo).

Hace exactamente lo que exige el skill:
  - mediana (no media)
  - exón-a-exón
  - exige >=3-5 exones CONSECUTIVOS para creer una deleción
  - declara el limite de deteccion
"""
import subprocess, statistics, sys

BAM = "/home/sw-cyborg/mva_hackathon/b_larga/bam/region_38_43M.bam"
# NM_001211 BUB1B, chr15, +, 23 exones
EXONS = [
 (40161068,40161255),(40165052,40165196),(40170061,40170121),(40170536,40170681),
 (40176476,40176673),(40183713,40183883),(40185164,40185379),(40185550,40185642),
 (40196544,40196774),(40199614,40199727),(40200243,40200359),(40200930,40200980),
 (40202404,40202465),(40202588,40202694),(40206183,40206458),(40208636,40208770),
 (40209634,40209775),(40210109,40210210),(40212498,40212648),(40213331,40213474),
 (40217495,40217667),(40218455,40218562),(40220563,40221122),
]

def median_depth(chrom, start, end):
    """Mediana de profundidad en la region (samtools depth -a incluye ceros)."""
    r = subprocess.run(
        ["samtools","depth","-a","-r",f"{chrom}:{start}-{end}",BAM],
        capture_output=True, text=True)
    deps = []
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) >= 3:
            deps.append(int(p[2]))
    if not deps:
        return None, 0, 0
    deps.sort()
    med = deps[len(deps)//2]
    return med, min(deps), max(deps)

print("=== CNV EXONICO BUB1B (NM_001211) — region_38_43M.bam (4 lanes, GRCh38) ===\n")
res = []
for i,(s,e) in enumerate(EXONS, 1):
    med, mn, mx = median_depth("chr15", s, e)
    res.append((i, s, e, e-s+1, med, mn, mx))
    print(f"exon {i:2d}  {s}-{e}  ({e-s+1:>4d} bp)  mediana={med}  min={mn}  max={mx}")

meds = [r[4] for r in res if r[4] is not None]
gmed = statistics.median(meds)
print(f"\nMEDIANA GLOBAL del gen: {gmed}")
print(f"(mediana de las medianas exonales)\n")

print("=== RATIO POR EXON (cobertura esperada diploide = 1.00) ===")
print(f"{'exon':>4} {'bp':>5} {'med':>5} {'ratio':>7}  {'flag'}")
for i,s,e,ln,med,mn,mx in res:
    ratio = med/gmed if med else 0
    flag = ""
    if ratio < 0.70: flag = "⚠️ <0.70 posible delecion het"
    elif ratio < 0.85: flag = "· <0.85 leve"
    elif ratio > 1.30: flag = "⚠️ >1.30 posible duplicacion"
    print(f"{i:>4} {ln:>5} {med:>5} {ratio:>7.3f}  {flag}")

print("\n=== RACHAS CONTIGUAS (>=3 exones) — el unico patron creible ===")
run = []
for i,s,e,ln,med,mn,mx in res:
    ratio = med/gmed if med else 0
    if ratio < 0.70:
        run.append((i, ratio))
    else:
        if len(run) >= 3:
            print(f"  ⚠️ RACHA de {len(run)} exones bajos: exones {run[0][0]}-{run[-1][0]}")
        run = []
if len(run) >= 3:
    print(f"  ⚠️ RACHA de {len(run)} exones bajos: exones {run[0][0]}-{run[-1][0]}")
else:
    print("  (sin rachas de >=3 exones contiguos bajo 0.70)")

minr = min(r[4]/gmed for r in res if r[4])
maxr = max(r[4]/gmed for r in res if r[4])
print(f"\nRango de ratios: {minr:.3f} — {maxr:.3f}")
print(f"\nLIMITE DE DETECCION: con mediana {gmed}x y exones de {min(r[3] for r in res)}-{max(r[3] for r in res)} bp,")
print("una delecion het de 1 exon daria ratio ~0.5 (caida del 50%) — DETECTABLE.")
print("Sin validacion ortogonal (MLPA/array-CGH).")
