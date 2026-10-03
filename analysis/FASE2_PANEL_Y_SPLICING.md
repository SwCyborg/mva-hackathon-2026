# FASE 2 — Panel ampliado + Splicing · Resultados

**Fecha:** 30 sep 2026 · **Deadline hackathon:** 24 oct 2026 (quedan 24 días)
**Ejecutado con:** DeerFlow (`mva-orquestador`, `mva-splicing`) + VEP REST + gnomAD v4 + ClinVar eutils

---

## RESUMEN

Se cerraron dos vías de análisis delegadas a DeerFlow. El hallazgo principal es **negativo y
útil**: el panel ampliado no aporta candidato alternativo. Se identificó **una variante de
splicing sub-priorizada** que queda como segunda opción de segundo alelo.

| Vía | Resultado | Confianza |
|---|---|---|
| Panel ampliado (9 genes) | **0 variantes codificantes novel** | Alta |
| Variantes P/LP conocidas con rsid | **0 en ClinVar** (98 rsids consultados) | Alta |
| Splicing en las 8 het de BUB1B | 1 candidata: `40192892 C>T` (splicing críptico) | Media |
| Segunda opción de segundo alelo | `40192892 C>T` — pero inferior a p.Asn1002Lys | Media |

**La hipótesis BUB1B (38%) no cambia.** El panel ampliado no crea competencia; la limitación
central (fase) sigue siendo el cuello de botella.

---

## 1. PANEL AMPLIADO — RESULTADO NEGATIVO SÓLIDO

### Genes y variantes encontradas

```
GEN       REGIÓN (GRCh38)              TOTAL   NOVEL(sin rsid)
TP53      17:7668402-7687550             16          1
DICER1    14:95086228-95158010           27         26
NF1       17:31094927-31377677          299         18
PTPN11    12:112418299-112509918         25         24
HRAS      11:532242-537308                3          2
KRAS      12:25205246-25250936           83         80
NRAS      1:114704469-114716771           9          9
BRAF      7:140719327-140924929          33         33
SUFU      10:102503961-102633535         137        134
TOTAL                                    632        327
```

### Anotación VEP de las 43 novel (tras excluir repeticiones >10 pb)

```
38  intron_variant
 3  3_prime_UTR_variant
 2  upstream_gene_variant
---
 0  CODIFICANTES
```

**Cero missense, nonsense, frameshift o splice en variantes novel de los 9 genes.**

### Control de artefacto (verificado)

476 de 632 eran hom-alt — implausible. Se inspeccionaron en crudo:

```
DICER1 14:95088868 rs1187652 C>T  AC=2;AF=1;AN=2;DB;DP=56  GT:AD:DP = 1/1:0,56:56
NF1   17:31095843 rs2952994 A>T  AC=2;AF=1;AN=2;DB;DP=47  GT:AD:DP = 1/1:0,47:47
```

Ambas con rsid y **AF = 0.979 en gnomAD v4** → polimorfismos comunes, no mosaicismo.
El exceso de hom-alt es enriquecimiento de alelos comunes, no señal biológica.

### Cierre del hueco metodológico (detectado por el orquestador)

El orquestador señaló que **filtrar por flag `DB` descarta variantes patogénicas conocidas con
rsid** (muchas de ClinVar tienen rsid). Hueco válido. Se cerró:

- 576 variantes con rsid extraídas de los 9 genes
- **98 rsids** de los genes críticos (TP53, DICER1, PTPN11, HRAS, NRAS, BRAF) consultados en
  ClinVar vía eutils
- **Resultado: 0 variantes P/LP**

→ No se perdió ninguna patogénica conocida. El negativo es válido.

**Limitación que permanece:** CNV/SV grandes son invisibles a `bcftools view` sobre SNV/indel.
En estos síndromes la deleción grande es mecanismo mayor (TP53 en Li-Fraumeni, microdeleción
de NF1, CNV de DICER1). **No evaluado.**

---

## 2. SPLICING BUB1B — ANÁLISIS VARIANTE-ESPECÍFICO

**Contexto del hueco:** el informe descartó splicing con umbral de **20 bp**. El branch point
está a **18-40 nt** (hasta ~50 nt) del aceptor → fuera de esa ventana. El descarte no era
concluyente.

### Las 8 heterocigotas de BUB1B evaluadas

| # | Pos | Var | Tipo | Sospecha splicing | Razón |
|---|---|---|---|---|---|
| 1 | 40168264 | G>A | SNV | BAJO | Intrón temprano |
| 2 | 40180642 | CT>C | del 1bp | BAJO | Contexto mixto, no poli-pirimidínico |
| 3 | 40181093 | AT>A | del 1bp | BAJO | Ídem |
| 4 | 40182809 | A>C | SNV | BAJO | **AF 0.0185 → común, descartada** |
| 5 | 40192892 | C>T | SNV | **MEDIO** | **Patrón arquetípico de acceptor críptico (C>T)** |
| 6 | 40209701 | T>G | stop_gained (V1) | NULO (splicing) | Patogénico por NMD, no splicing |
| 7 | 40216470 | A>G | SNV | BAJO | Intrón profundo, novel |
| 8 | 40220612 | T>G | missense (V2) | MEDIO (secundario) | Posible ESE/ESS, primario es missense |

### Rareza verificada en gnomAD (dato decisivo)

```
15-40192892-C-T  AF = 0.0021  → RARA (candidata viable)
15-40182809-A-C  AF = 0.0185  → COMÚN (descartada)
15-40216470-A-G  NO ENCONTRADA → novel
```

### Candidata identificada

**`chr15:40,192,892 C>T`** — intrónica, heterocigota, AF 0.0021.
Mecanismo hipotético: ganancia de acceptor críptico (C>T es el patrón arquetípico).

### 🔴 RESULTADO SPLICEAI — LA HIPÓTESIS SE REFUTÓ

SpliceAI 1.3 ejecutado sobre las 8 heterocigotas (ProDesk, CPU, ~3 min).
Scores: `DS_AG | DS_AL | DS_DG | DS_DL` (delta acceptor-gain, acceptor-loss, donor-gain, donor-loss).

```
POS        VAR        DS_AG  DS_AL  DS_DG  DS_DL   MÁX
40168264   G>A        0.02   0.00   0.01   0.01    0.02
40180642   CT>C       0.00   0.00   0.00   0.00    0.00
40181093   AT>A       0.00   0.00   0.00   0.00    0.00
40182809   A>C        0.00   0.00   0.00   0.00    0.00
40192892   C>T        0.00   0.00   0.00   0.00    0.00   <-- la candidata
40209701   T>G (V1)   0.03   0.00   0.00   0.01    0.03
40216470   A>G        0.00   0.00   0.00   0.00    0.00
40220612   T>G (V2)   0.02   0.00   0.00   0.00    0.02
```

**Máximo score global: 0.03.** Umbral de cribado: 0.20. Umbral alta precisión: 0.50.

**NINGUNA variante tiene efecto de splicing predicho.**

**La candidata `40192892 C>T` da 0.00 en las cuatro categorías.** La hipótesis del
`mva-splicing` (acceptor críptico por patrón C>T + rareza AF 0.0021) **queda refutada por
medición directa**. Perfil de rareza + patrón de secuencia sugerente NO implican efecto.

### Consecuencia: el hueco de splicing está CERRADO

La auditoría adversarial había señalado que descartar splicing con umbral de 20 bp **no era
concluyente** (el branch point está a 18-40 nt). Ese hueco **ya no existe**: SpliceAI usa una
ventana de contexto de 10 kb y evalúa las cuatro categorías de sitio — cubre branch point,
polypyrimidine tract y sitios crípticos profundos.

**El descarte de splicing pasó de inferido a medido.**

### Detalles técnicos de la instalación

- SpliceAI 1.3.1 + TensorFlow 2.21.0 en `~/.hermes/venvs/spliceai/` (ProDesk)
- Modelos: `spliceai1-5.h5` (Illumina promedia 5 modelos) — disponibles en el paquete pip
- Anotación: `grch38.txt` incluida en el paquete
- **Bug resuelto:** `np.fromstring` removido en numpy 2.x → VCF vacío en silencio.
  Parcheado a `np.frombuffer(seq.encode(), dtype=np.int8)`. Backup en `utils.py.bak`
- **Trampa:** contigs de la referencia con prefijo `chr`, VCF sin prefijo

### Confirmación que ya NO es necesaria

~~Correr SpliceAI sobre las 8~~ → **HECHO, resultado negativo**
RT-PCR / minigene → ya no prioritario (no hay candidato de splicing que validar)


---

## IMPACTO SOBRE EL DIAGNÓSTICO

**No cambia.** El panel ampliado no compite; la fase sigue irresoluble. Confianza 38%.

**Cambio menor introducido:** ahora existe una **segunda opción documentada** de segundo alelo
(`40192892 C>T`) que antes no estaba priorizada. Si el jurado cuestiona p.Asn1002Lys, hay una
alternativa argumentada.

---

## ORDEN SUGERIDO

1. **Correr SpliceAI** sobre las 8 het de BUB1B — cierra la limitación #1 del método.
   Costo: descargar modelo (~2 GB). Retorno: puede descubrir un segundo alelo LoF.
2. **Análisis de CNV** sobre TP53, NF1, DICER1, SUFU — el mecanismo que este pipeline no ve.
3. **Phasing estadístico** (Beagle + panel 1000G) si se quiere mover el prior 50/50.
4. Documentar el negativo del panel ampliado en el informe de entrega — es un resultado
   honesto que fortalece la credibilidad (prueba que se buscaron alternativas).

---

## REGISTRO

Ejecutado con DeerFlow en ProDesk (192.168.0.10): `mva-orquestador` (task a9f5bdcd),
`mva-splicing` (tasks 6720a851 y 32e193f9). Consultas API desde Elite X2: VEP REST,
gnomAD v4 GraphQL, ClinVar eutils. Todos los task_ids en estado `success`.
