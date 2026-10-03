# SwCyborg — MVA Hackathon 2026, Track 1 (Variant Prediction)

Código reproducible, submission y reporte para el **Rare Disease, Real Kid: MVA Hackathon 2026**.

**Equipo:** SwCyborg · **Proband:** PROBAND01

---

## Resultado

**Candidato causal: BUB1B compound-het** (MVA1, OMIM 257300)

| # | Variante (GRCh38) | Consecuencia | Evidencia |
|---|---|---|---|
| V1 | `chr15:40,209,701 T>G` | stop_gained p.Leu737Ter | ClinVar Pathogenic/LP (VCV000533901), trait MVA1 exacto |
| V2 | `chr15:40,220,612 T>G` | missense p.Asn1002Lys | SIFT 0.01, PolyPhen 0.997, phyloP 4.80, novel |

**Confianza del genotipo como diagnóstico: 38%** = `0.95 × 0.80 × 0.50`.
El factor 0.50 es la **fase cis/trans sin resolver** — limitación dura del short-read sin trío.

---

## Por qué este repo es reproducible

Todo el pipeline está en `scripts/` y cada número del reporte se puede regenerar desde el VCF crudo.

```bash
# 1. Verificar las variantes candidatas + extraer los 7 genes MVA
python3 scripts/pipeline_analysis.py --vcf <WGS_*.vcf.gz> --out analysis/

# 2. Validar el submission contra la lógica del evaluador oficial
python3 scripts/validar_submission.py submission/Swcyborg_track1_submission.csv

# 3. Cierre del CNV exónico de BUB1B (requiere BAM alineado GRCh38)
python3 scripts/cnv_bub1b_exones.py
```

---

## Estructura

```
scripts/
  pipeline_analysis.py       Pipeline completo: verificación de variantes + control negativo
  validar_submission.py      Replica la lógica de evaluation.py del evaluador oficial
  cnv_bub1b_exones.py        Llamado por exón del CNV de BUB1B (23 exones)
submission/
  Swcyborg_track1_submission.csv   CSV final (5 filas, epcr descendente)
report/
  Swcyborg_track1_report.md        Reporte de métodos (requerido para jueces)
analysis/
  (generado)                  Salidas JSON del pipeline
```

---

## Hallazgos metodológicos clave

### 1. Los 5 mecanismos alternativos fueron evaluados y cerrados

| Mecanismo | Veredicto | Evidencia |
|---|---|---|
| CNV exónico BUB1B | **Negativo con poder** | Mediana 54×, ratios 0.722-1.463, 0 bases DP<10, límite ~1 exón |
| Splicing BUB1B | Refutado | SpliceAI Δscore máx 0.03 vs umbral 0.20 |
| Mosaicismo puntual | Refutado | Ver `report/` (refutación documentada) |
| Otros 6 genes MVA | Excluidos | Sin bloques contiguos de profundidad reducida (DP VCF 39-52×) |
| Fase cis/trans | **Irresoluble (medido)** | 0 lecturas unen V1↔V2, distancia 10.9 kb vs inserto 50 bp |

### 2. La trampa del CNV con referencia acotada

Un intento previo alineó contra una referencia acotada de 863 KB (9 genes). Los BAM resultantes
son **inservibles para CNV** — todo read con homología parcial mapea ahí obligado, inflando la
cobertura artificialmente (TP53 mediana 32× con `max/mediana` = 5.818×, biológicamente imposible).
El CNV válido requiere alineación contra GRCh38 completo. Documentado en `scripts/cnv_bub1b_exones.py`.

### 3. Nomenclatura de cromosoma

El VCF del dataset usa contigs **sin prefijo `chr`** (`##contig=<ID=15,...>`), mientras el template
de submission oficial usa `chr15`. El evaluador compara `chrom` como string sin normalizar →
la ambigüedad es binaria. El submission cubre ambas nomenclaturas.

---

## Declaración de uso de IA

- **DeepSeek V4 Pro** (DeepSeek via OpenRouter API): orquestado por el runtime Hermes Agent
  ("Iris"). Términos comerciales; **sin entrenamiento sobre contenido del cliente**.
- **Anthropic Claude** (Claude Code CLI): usado para **revisión adversarial de la estrategia de
  submission** y depuración del pipeline — **no** para generar anotaciones de variantes.
- **Todas las anotaciones** (VEP, ClinVar, gnomAD, UCSC, SpliceAI) fueron computadas con
  **herramientas y APIs públicas**, no generadas por un LLM.

---

## Datos

Los datos genómicos del probando **NO están incluidos** en este repo (acceso controlado, CC BY 4.0
en las salidas). Descargables desde `SageBio/mva-hackathon-2026-data` previa aceptación de los
términos del hackathon.

## Licencia

Salidas publicadas bajo **CC-BY 4.0** conforme a las reglas del hackathon. Sin intento ni intento
de re-identificación. No se contactó al probando ni a su familia.
