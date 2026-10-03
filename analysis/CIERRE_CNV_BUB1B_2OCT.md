# CIERRE HUECO #5 — CNV exónico de BUB1B

**Fecha:** 2 de octubre de 2026
**Estado:** ✅ CERRADO — negativo con poder estadístico
**Confianza global MVA:** 38% (sin cambio; la fase sigue siendo el limitante)

---

## 1. Resultado

**No se detecta deleción ni duplicación en ninguno de los 23 exones de BUB1B.**

```
BUB1B · NM_001211 · chr15:40.161.068-40.221.122
BAM: b_larga/bam/region_38_43M.bam  (4 lanes, GRCh38 completo)
Mediana global: 54× (mediana de las medianas exonales)
Rango de ratios exonales: 0.722 — 1.463
Rachas contiguas ≥3 exones bajo 0.70: NINGUNA
```

### Detalle por exón

| Exón | bp | mediana | ratio | observación |
|---|---|---|---|---|
| 1 | 188 | 44 | 0.815 | ruido |
| 2 | 145 | 42 | 0.778 | ruido |
| 3 | 61 | 56 | 1.037 | — |
| 4 | 146 | 45 | 0.833 | ruido |
| 5 | 198 | 62 | 1.148 | — |
| 6 | 171 | 46 | 0.852 | — |
| 7 | 216 | 60 | 1.111 | — |
| 8 | 93 | 46 | 0.852 | — |
| 9 | 231 | 54 | 1.000 | — |
| 10 | 114 | 55 | 1.019 | — |
| 11 | 117 | 53 | 0.981 | — |
| 12 | 51 | 68 | 1.259 | ventana chica (Poisson) |
| 13 | 62 | 79 | 1.463 | ventana chica (Poisson) |
| 14 | 107 | 43 | 0.796 | ruido |
| 15 | 276 | 53 | 0.981 | — |
| 16 | 135 | 50 | 0.926 | — |
| 17 | 142 | 55 | 1.019 | — |
| 18 | 102 | 67 | 1.241 | ventana chica |
| 19 | 151 | 63 | 1.167 | — |
| 20 | 144 | 58 | 1.074 | — |
| 21 | 173 | 39 | 0.722 | ruido |
| 22 | 108 | 66 | 1.222 | ventana chica |
| 23 | 560 | 51 | 0.944 | — |

---

## 2. Por qué el negativo SOSTIENE (y no es un negativo sin poder)

Regla del skill: un negativo solo vale si el test TIENE PODER. Los cuatro
indicadores que lo confirman:

| Indicador | Valor | Umbral | Veredicto |
|---|---|---|---|
| Mediana de cobertura | **54×** | ≥20-30× (industria) | ✅ sobra |
| Bases DP=0 | **0** (0,00%) | ~0% | ✅ sin huecos |
| Bases DP<10 | **0** (0,00%) | bajo | ✅ toda la región con poder |
| `DP_max / mediana` | **1,76** | <4 (sano) | ✅ distribución sana |
| Exones evaluados | **23/23** | 100% | ✅ cobertura completa del gen |

**Poder de detección:** con mediana 54× y exones de 51-560 bp, una **deleción
heterocigota de un solo exón** daría ratio ~0.5 (caída del 50%). El ruido
observado (0.722-1.463) está muy lejos de ese valor. **El test puede excluir
una deleción de un exón.**

### Frase de cierre obligatoria

> *Negativo a 54× de cobertura (mediana), límite de detección ~1 exón, sin
> validación ortogonal (MLPA/array-CGH).*

---

## 3. La señal que NO se explica (declarada, no escondida)

Los exones **12 (1.259), 13 (1.463), 18 (1.241) y 22 (1.222)** están por
encima de 1.20. Son los **más pequeños** (51, 62, 102 y 108 bp) → mayor
varianza de Poisson por conteo bajo de bases.

**No constituye duplicación:**
1. **Sin racha contigua.** Una duplicación real afectaría un **bloque** de
   exones consecutivos. Aquí los elevados están dispersos (12-13, luego 18,
   luego 22) y el mínimo creíble del skill es **3-5 exones contiguos**.
2. **Va en la dirección segura.** El sesgo es *hacia arriba*, no hacia abajo.
   Un artefacto que infla no puede esconder una deleción (que resta señal).
3. **Explicable por diseño.** Ventanas cortas → menos bases → más varianza.
   Con 51 bp, la diferencia entre mediana 68 y 54 son ~14 reads: ruido normal.

---

## 4. El artefacto que se descartó por el camino

**Los BAMs del panel ampliado (`cnv_panel/aligned/*.panel.bam`) NO sirven para
CNV.** Se alinearon contra la referencia acotada `panel_genes.fa` (863 KB, 9
contigs con nombre de gen), no contra GRCh38.

Medición que lo demuestra (TP53, L001):

```
n=19.149 pb · p10=12 p25=16 MEDIANA=32 p75=3.459 p90=15.095 p99=99.909 max=186.137
bases DP>100: 8.768 (45,8%)
```

`p75/mediana = 108×`. Casi la mitad de las bases en la cola. **Cobertura media
calculada: 11.548×** — biológicamente imposible para un WGS.

**Causa (documentada en el skill):** al ofrecer solo 863 KB como referencia,
todo read con homología parcial mapea ahí **obligado** — no tiene alternativa.

| BAM | Referencia | Mediana | max/mediana | DP=0 | Veredicto |
|---|---|---|---|---|---|
| `cnv_panel/aligned/L001.panel.bam` | panel 863 KB | 32 | **5.818×** | 0 | ❌ artefacto |
| **`b_larga/bam/region_38_43M.bam`** | **GRCh38** | **54** | **1,76×** | **0** | ✅ **sano** |

**El panel ampliado de 9 genes NO puede cerrar ningún hueco CNV.** Su único
uso válido es comparar genes entre sí (y ahí ya rindió: HRAS es el de menor
target, no el de menor cobertura).

---

## 5. Método reproducible

```bash
# BAM SANO (GRCh38, 4 lanes)
BAM=~/mva_hackathon/b_larga/bam/region_38_43M.bam

# Test de sanidad ANTES de interpretar
samtools depth -a -r chr15:40155984-40226137 $BAM | cut -f3 | sort -n > /tmp/dp.txt
awk '{a[NR]=$1} END{printf "n=%d p10=%d p25=%d MEDIANA=%d p75=%d p90=%d max=%d\n",
     NR, a[int(NR*.1)], a[int(NR*.25)], a[int(NR*.5)], a[int(NR*.75)],
     a[int(NR*.9)], a[NR]}' /tmp/dp.txt
awk '$1<10' /tmp/dp.txt | wc -l     # bases sin poder  -> debe ser ~0
awk '$1>100' /tmp/dp.txt | wc -l    # cola             -> debe ser ~0

# Ratio por exón (script completo)
python3 /tmp/cnv_bub1b_exones.py
```

**Reglas aplicadas del skill:**
- Mediana, **nunca media** (la media miente con distribuciones de cola).
- Ordenar a archivo y luego indexar percentiles (no ordenar dentro del awk que agrega).
- `DP_max/mediana` < 4 y `%DP=0` ≈ 0 como criterio de sanidad.
- Exigir **≥3-5 exones contiguos** para creer una deleción.
- Declarar el límite de detección y la falta de validación ortogonal.

---

## 6. Impacto en el entregable del hackathon

**Huecos del caso y su estado (deadline 24 oct 2026):**

| Hueco | Estado | Evidencia |
|---|---|---|
| Fase cis/trans | ✅ Cerrado — irresoluble **medido** | 0 lecturas V1↔V2, distancia 10,9 kb vs inserto 50 pb |
| **CNV exónico BUB1B** | ✅ **CERRADO (este documento)** | 54× mediana, ratios 0.72-1.46, sin rachas |
| Panel ampliado (9 genes) | ✅ Cerrado | 0 codificantes novel; 98 rsids ClinVar → 0 P/LP |
| Splicing BUB1B | ✅ Cerrado — refutado por medición | SpliceAI Δscore máx 0.03 vs umbral 0.20 |
| Mosaicismo | ✅ Cerrado — refutado | (v. REFUTACION_MOSAICO.md) |

**Candidato:** BUB1B compound-het, **confianza 38%** = `0,95 × 0,80 × 0,50`.
El factor 0,50 es la **fase sin resolver** — el único hueco que queda abierto
y que ningún cómputo sobre los datos actuales puede cerrar (limitación dura
del short-read).

**Acción pendiente en el entregable:** el CSV (`Swcyborg_track1_submission.csv`)
y el reporte (`Swcyborg_track1_report.md`) **no mencionan este cierre**. Deben
actualizarse para reflejar que el negativo de CNV ahora tiene poder declarado.

⚠️ **Recordar la regla del 30 sep:** cada vez que una auditoría o un cierre
cambien una conclusión numérica, **regenerar el entregable en el mismo
movimiento**. Un artefacto generado una vez y consumido después guarda un
estado que venció.

---

## 7. Registro

- **Comando de medición:** `samtools depth -a -r chr15:40155984-40226137 region_38_43M.bam`
- **Script de exones:** `/tmp/cnv_bub1b_exones.py`
- **BAM:** `~/mva_hackathon/b_larga/bam/region_38_43M.bam` (101 MB, mtime 30 sep 08:43)
- **Exones:** `ref/refFlat.txt` → `NM_001211` (23 exones)
- **Predecesor refutado:** `VERIF_CNV_BUB1B.md` (26 sep — el negativo viejo sin poder,
  sobre 1 lane y bins agregados)
