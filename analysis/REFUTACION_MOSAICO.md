# REFUTACIÓN — Tu hallazgo de aneuploidía mosaica NO resiste el control

## Qué propusiste

Priorizaste formalizar el screening de mosaicismo (chr16/20/21/22) con MoChA como
evidencia #1 para el jurado. Dato original:

```
chr16 mean_DP=58.83 (+37%)
chr20 mean_DP=49.40 (+18%)
chr21 mean_DP=56.85 (+35%)
chr22 mean_DP=51.51 (+22%)
resto autosomas: 41-43x
chrX: 21.79x (~0.5x autosomal)
```

Con el caveat correcto: "exploratorio, no concluyente... sesgo de GC, mapeabilidad,
selección de sitios del caller... chr16/21/22 son acrocéntricos con repeticiones rDNA
que pueden inflar profundidad por mismapeo."

**Tenías razón en poner el caveat. Lo que faltaba era dimensionarlo.** Lo corrí como test.

## Reproduje tus números exactos (1,5M variantes, no muestra)

```
chr1   43,75   n=383.428
chr2   43,58   n=388.761
chr15  43,99   n=140.405
chrX   22,53   n=117.225   ← 0,52x autosómico
chr16  60,30   n=143.762   ← +38%
chr20  51,13   n=119.444   ← +17%
chr21  58,53   n=86.333    ← +34%
chr22  52,73   n=83.380    ← +21%
```

Tus valores se confirman. El promedio cromosómico es correcto. **Pero el promedio es
la medida equivocada para distinguir dosis de mapeabilidad.**

## Test 1 — DP por región (p-arm vs q-arm)

Una trisomía eleva la profundidad UNIFORMEMENTE en todo el cromosoma.
El exceso localizado en una clase de secuencia es la firma del artefacto.

```
chr16   p-arm(0-15Mb):  44,02  │  q-arm(>15Mb):  65,40
chr20   p-arm(0-15Mb):  43,33  │  q-arm(>15Mb):  53,19
chr21   p-arm(0-15Mb):  95,05  │  q-arm(>15Mb):  43,19   ← invertido
chr22   p-arm(0-15Mb):  93,77  │  q-arm(>15Mb):  43,39   ← invertido
chr2    p-arm(0-15Mb):  42,74  │  q-arm(>15Mb):  43,64   ← control plano
```

El exceso NO está distribuido. Está confinado a ~15 Mb de un solo brazo, y en
direcciones opuestas. chr21/22 → brazo corto (rDNA). chr16/20 → brazo largo
(heterocromatina pericentromérica; chr16 tiene el mayor bloque del genoma, 16q11.2).

**Un promedio sobre 0-250 Mb esconde que todo el exceso viene de ~6% del cromosoma.**
Eso es exactamente el sesgo de selección de sitio que mencionaste — pero no es un
desperfecto menor, es el 100% de la señal.

## Test 2 — Distribución de BAF (mecanismo, no correlación)

Este es el test que puede confirmar o matar el hallazgo de forma mecanística.
Trisomía desplaza la BAF de hets hacia 0,33/0,67. Conté la fracción de hets con
BAF en 0,45-0,55 (diploide esperado: ~95%+).

**Control positivo primero — para validar que el método SÍ detecta cambio de dosis:**
```
chrX   n=68.118   BAF~0,5 = 3,3%    ← hemicigoto varón: casi sin hets reales
chr1   n=348.849  BAF~0,5 = 27,1%
chr2   n=364.044  BAF~0,5 = 27,6%
```
chrX a 3,3% cuando el valor diploide es ~27%: **el método detecta pérdida de dosis.**

**Ahora los sospechosos:**
```
chr16  26,9%   ← +38% DP
chr20  26,1%   ← +17% DP
chr15  26,8%   ← control
chr1   27,1%   ← control
chr2   27,6%   ← control
```

Los 4 sospechosos tienen la MISMA distribución de BAF que los controles sanos.
Si chr16 fuera trisómico, su fracción BAF~0,5 caería. No cae.

## Test 3 — BAF por brazo: el mecanismo queda expuesto

```
chr21  p-arm: BAF~0,5 = 11,2%   ← colapsada
chr21  q-arm: BAF~0,5 = 26,7%   ← normal

chr22  p-arm: BAF~0,5 = 10,4%   ← colapsada
chr22  q-arm: BAF~0,5 = 27,1%   ← normal

chr16  p-arm: 26,4%  │  chr16 q-arm: 27,1%   (plano)
chr20  p-arm: 30,9%  │  chr20 q-arm: 24,8%   (plano)
```

El brazo corto de chr21/22 tiene DP 95x y BAF colapsada al mismo tiempo. Esa es la
firma inequívoca de reads mal mapeados apilados sobre repeticiones: reads de muchas
regiones del genoma aterrizan en el rDNA, suben la profundidad, y generan "hets" con
AD artificialmente balanceados que no son heterocigotos reales.

chr16/20 muestran el mismo efecto en grado menor y distribuido (heterocromatina).

## Veredicto

**REFUTADO.** No hay aneuploidía mosaica detectable. El hallazgo era la suma de dos
artefactos de mapeabilidad:

1. rDNA en brazos cortos acrocéntricos (chr21/22) — exceso localizado 2,2x
2. Heterocromatina pericentromérica (chr16/20) — exceso difuso moderado

**Consecuencia para tu recomendación:** la prioridad #1 que propusiste para el jurado
("formalizar con MoChA") queda anulada. MoChA usa BAF+LRR sobre VCF — habría formalizado
el mismo artefacto con más aparato estadístico y más credibilidad aparente. Eso es peor
que no hacerlo: un artefacto con sello de herramienta validada es más difícil de retractar.

## Lo que TUS hallazgos sí aportaron (y son valiosos)

1. **El fasing por lecturas es IMPOSIBLE, no difícil.** Yo reporté "VCF sin PS/PID"
   (correcto pero subdimensionado). Tú mostraste que con solo 3 hets en 25 kb y ningún
   PID compartido, WhatsHap no tiene cadena que armar. Eso mató las Opciones A/B/C y
   me ahorró 20-40 h de cómputo. **Ese hallazgo es el más valioso del día.**

2. **El filtrado por k-mero con BBDuk** (2-4 h, <2 GB) en vez de alinear 85 GB.
   Correcto: mi estimado de 20-40 h y la inviabilidad de disco eran acertados, pero
   la vía barata no la había considerado.

3. **Priorizar transparencia sobre la fase** como argumento metodológico ante un
   jurado cuantitativo. Correcto y lo adopto.

## Nota de método

Corrí el BAF dos veces: la primera falló por un bug de quoting mío (`printf` con
`<0.35` interpretado como redirección) y por usar `asort` (gawk-only). Lo corregí
y re-corrí. Reporto solo la corrida válida.

Todo sobre ProDesk (100.95.59.11), solo lectura sobre el VCF y el BAM existente.
No modifiqué datos.
