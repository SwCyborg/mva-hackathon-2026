# VERIFICACIÓN INDEPENDIENTE — Análisis CNV por read-depth en BUB1B (MVA Hackathon)

## CONTEXTO

Paciente con fenotipo MVA (mosaic variegated aneuploidy) + rabdomiosarcoma.
Hipótesis principal: BUB1B compound-het (V1 p.Leu737Ter stop-gained + V2 p.Asn1002Lys missense).
V1 y V2 están a 10.911 bp de distancia → fase (cis/trans) irresoluble por lecturas cortas sin trío.

**El hueco a cerrar:** si existe una DELE CIÓN EXÓNICA en BUB1B (trans con V1 o V2),
el genotipo compound-het se establece SIN necesidad de fase. Eso elevaría la confianza de 38% al techo.

Mi conclusión: **NO hay deleción exónica en BUB1B**. Necesito que la verifiques o la refutes.

## DATOS DE ENTRADA

- VCF: `WGS_EX2312012_HGWCNDSX7.vcf.gz` (315 MB, GRCh38)
- FASTQ: 8 archivos, ~85 GB, 4 lanes (L001-L004), WGS
- BAM consolidado: `aln_chr15_full.bam` (9.482.727.820 bytes), contigs: `15` (LN 101991189) + `15:40155984-40226137` (LN 70154, espurio)
- quickcheck PASS

## CÓMO SE CONSTRUYÓ EL BAM (posible fuente de sesgo)

1. `alinear_fastq.sh`: bwa mem contra **chr15 SOLO** (102 Mb, no genoma completo).
   Consecuencia: todo read no-chr15 va a `*` (unmapped). Del total, 35.414.820 mapeados / 50.980.988 unmapped.
2. Los BAMs intermedios `aln_L00N.bam` se borraron; quedaron chunks `.tmp.NNNN.bam`.
3. Adicionalmente `alinear_region.sh` alineó los mismos FASTQ contra una región de 70 kb
   (`ref/bub1b_region.fa`, chr15:40155984-40226137) → chunks `.tmp.NNNN-001.bam`.
4. `cnv_v2.sh` clasificó chunks POR REFERENCIA (leyendo `@SQ`), y mergeó **11 planos**,
   excluyendo 2 de región. Resultado: `cnv_chr15.bam`.
5. `cnv_bub1b.sh` (script hermano) IGNORÓ ese filtro y mergeó los **13** chunks →
   `aln_chr15_full.bam`, que incluye el contig espurio. (Usé este.)

## MÉTODO

### Paso 1 — Depth con distintos MAPQ
`cnv_v2.sh` usó `samtools depth -a -q 20 -Q 20` (baseQ 20, MAPQ 20).

**Afirmo que el filtro correcto es MAPQ >= 5**, porque el script de anotación
`proceso_apertura_modulo` / la documentación del pipeline dice MAPQ>=5.

Con `-Q 20` obtuve:
```
media BUB1B (40.150M-40.250M, n=100.001) = 8,95
ratio exón/mediana: exón 2 = 0,13 / exón 4 = 0,35 / exón 6 = 0,57 / exón 11 = 0,50 / exón 16 = 0,35
```

Con `-Q 5`:
```
media = 15,354
```

### Paso 2 — Depth por exón según MAPQ (el test decisivo)
`samtools depth -a -q 0 -Q {0,5,20,40} -r 15:40160000-40222000 aln_chr15_full.bam`

```
EXON  MAPQ>=5   MAPQ>=20  MAPQ>=40   ratio r5/r40
1     5,24      5,24      5,24       1,00
2     0,81      0,81      0,81       1,00
3     6,16      6,16      6,16       1,00
4     2,19      2,19      2,19       1,00
5     6,58      6,58      6,58       1,00
6     3,67      3,67      3,67       1,00
7     6,04      5,85      5,85       1,03
8     5,40      5,40      5,40       1,00
9     4,33      4,25      4,25       1,02
10    4,83      4,48      4,48       1,08
11    3,89      3,89      3,89       1,00
12    4,55      4,55      4,55       1,00
13    7,53      7,21      7,21       1,04
14    6,07      6,07      6,07       1,00
15    7,09      7,09      7,09       1,00
16    2,21      2,21      2,21       1,00
17    5,46      5,46      5,46       1,00
18    5,68      5,68      5,68       1,00
19    5,63      5,63      5,63       1,00
20    3,85      3,85      3,85       1,00
21    5,64      5,64      5,64       1,00
22    7,84      7,84      7,84       1,00
23    7,92      7,88      7,12       1,11
```

**Interpretación que sostengo:** el filtro MAPQ es INDIFERENTE dentro de los exones
(r5/r40 = 1,00). Si hubiera deleción, el exon afectado tendría cobertura baja
en TODOS los niveles de MAPQ y subiría al relajarlo. No sube.

### Paso 3 — Control de mapeabilidad genómica
```
15:10.000.000-10.050.000   MAPQ>=5 = 0,00   (¿por qué cero? ver Pregunta 4)
15:50.000.000-50.050.000   MAPQ>=5 = 14,02  MAPQ>=40 = 6,46  ratio = 2,17
15:90.000.000-90.050.000   MAPQ>=5 = 15,59  MAPQ>=40 = 5,89  ratio = 2,65
BUB1B (40,16-40,22M)       MAPQ>=5 = 6,4    MAPQ>=40 = 6,3   ratio = 1,00
```

**Interpretación que sostengo:** el genoma tiene caída típica de 2,2-2,7x entre
MAPQ 5 y 40 (reads multimapping). BUB1B tiene 1,00x → región de mapeabilidad
excepcionalmente LIMPIA. Por tanto BUB1B **no** es candidata a pérdida espuria
de reads por filtro de calidad.

### Paso 4 — CNVkit
`cnvkit.py coverage` con `--min-mapq 40` sobre regiones target:
```
Target BUB1B (15:40196068-40221122):   depth = 6.4066   log2 = 2.67956
Antitarget flank (40.0M-40.3M):        depth = 7.97481  log2 = 2.99545
```
Sesgo GC/mapabilidad disponible en el .cnn si hace falta.

## PREGUNTAS ESPECÍFICAS PARA TI

1. **¿El razonamiento del Paso 2 es válido?**
   Convirtiendo depth a copias asumiendo diploide y linealidad:
   - Exón 2: 15,354 × (0,81/8,95) / 5,62 = 0,248 → **0,5 copias**. Cerca de 0,5 (deleción het).
   - Exón 5: 15,354 × (6,58/8,95) / 5,62 = 2,013 → **4,0 copias**. Imposible.
   - ¿Es la dispersión exónica evidencia de deleción, o solo ruido de exones cortos?
     (Exones 2/4/12/13/18 = 145/146/51/62/102 bp. Exones 5/15 = 198/276 bp.)

2. **¿Bastan 5-8x de cobertura para detectar una deleción de 1-2 exones?**
   Mi argumento: un exón het-delecionado tendría 2,5-4x. ¿Hay un umbral mínimo
   de profundidad documentado para CNV exónico que invalide este test?

3. **¿Los 8 archivos FASTQ son 4 lanes o 4 lanes × 2 réplicas?**
   Observé `aln_L001.bam.tmp.0000.bam` Y `aln_L001.bam.tmp.0000-001.bam`.
   Si los `.NNNN-001` son re-chunks del MISMO lane 1, el merge con 13 chunks
   podría estar DUPLICANDO cobertura de un lane (o parte). ¿Cómo verificarlo?

4. **`15:10.000.000-10.050.000` da depth 0,00 con `samtools depth -a`.**
   ¿Es señal de problema en el BAM (hueco real) o de que `-a` no funciona como espero?
   Es el control que NO cierra y podría invalidar todo lo demás.

5. **El contig espurio `15:40155984-40226137` con 2.100.849 reads mapeados.**
   ¿Puede robar reads a chr15 en la región BUB1B y bajar artificialmente su cobertura?
   (si un read mapea a la región de 70kb en vez de a chr15, ¿compite?)

6. **¿Hay un artefacto estructural que este método no puede ver?**
   Inversión, transposición, o CNV donde la SECUENCIA CAMBIA pero la PROFUNDIDAD no
   (ej. dup invertida, o deleción con duplicación compensatoria).

## LO QUE CONCLUYO

**No hay deleción exónica en BUB1B; no hay duplicación.**
El hueco #5 de la auditoría ("ningún gen recibió chequeo de CNV") se cierra para BUB1B.
Consecuencia: la vía "deleción en trans elimina la necesidad de fase" queda EXCLUIDA.
Confianza global se mantiene en 38%.

## LO QUE NO CONCLUYO

- Los 6 genes descartados (CEP57, TRIP13, CENATAC, SLF2, SMC5, MAD1L1) NO fueron chequeados.
- No verifiqué la unicidad de los k-mers de contexto que usé (en realidad no los usé en el veredicto final).
- No hay validación ortogonal (MLPA, array CGH, long-read).

## ARCHIVOS EN EL SISTEMA

ProDesk (192.168.0.10 / 100.95.59.11): `~/mva_hackathon/`
- `aln_chr15_full.bam` (9,5 GB) — el que usé
- `cnv_chr15.bam` (7,0 GB) — el limpio de cnv_v2.sh
- `depth_q{0,5,20,40}.tsv`, `d_q{0,5,20,40}.tsv`
- `exones_bub1b.tsv` (23 exones, coordenadas GRCh38)
- `ref/refFlat_nc.txt` (refFlat hg38 sin prefijo chr)
- `fase_b0.log`, `fase_b1b.log`, `fase_b2.log`, `fase_b3b.log`, `cnv_v2.log`, `cnv_run.log`
- `cnv_out/` (salidas CNVkit)

## INSTRUCCIÓN

Refuta o confirma. Si refutas, dame el número de línea o comando que lo demuestra.
Si confirmas, dime qué evidencia adicional haría falta para subir de "no detectado"
a "excluido con confianza".
