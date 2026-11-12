# Doerl_2025

Companion analyses for:

> **Maturation of the glymphatic system confers innate resistance of the brain to Zika virus infection**
> Doerl JG, Chaves EKC, dos Santos VAB, Marcelino BLM, Coelho DM, Miyakoshi LM, Araujo JMG,
> Jeronimo SMB, Sequerra EB.
> bioRxiv 2025.10.26.684625 — [preprint](https://www.biorxiv.org/content/10.1101/2025.10.26.684625v1)
> · [doi:10.1101/2025.10.26.684625](https://doi.org/10.1101/2025.10.26.684625)
> · [preLights highlight](https://prelights.biologists.com/highlights/maturation-of-the-glymphatic-system-confers-innate-resistance-of-the-brain-to-zika-virus-infection/)

The preprint shows that a maturing glymphatic system acts as an innate barrier to ZIKV. Injected
before maturation, the virus spreads widely, deforms vessels and disrupts the Aqp4⁺ astrocyte
network; after maturation the same challenge causes neither disease nor persistent infection.
Blocking vasopressin receptor 1a, which reduces perivascular fluid dynamics, reopens the window of
susceptibility.

**📊 Reports: <https://mcblab.github.io/Doerl_2025/>**

## Analyses

| Folder | Dataset | What it asks |
|---|---|---|
| [`macaque/`](macaque/) | GEO [GSE342497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE342497) — infant rhesus macaque scRNA-seq, ZIKV vs control, 3 brain regions ([Edara et al., *JCI Insight* 2026](https://doi.org/10.1172/jci.insight.203088)) | What happens to AQP4 and the rest of the perivascular anchoring complex in a brain infected *after* the equivalent maturation window? [→ report](https://mcblab.github.io/Doerl_2025/macaque/) |

Headline from `macaque/`: AQP4 transcript is unchanged and there is no leukocyte influx — consistent
with a protected, mature brain — but the machinery that holds AQP4 at the vessel (`DAG1`, `KCNJ10`,
`MLC1`, `UTRN`, `LARGE1`, and the basement-membrane ligands `NID1`, `LAMB2`, `LAMC1`, `HSPG2`,
`AGRN`) is coordinately lower in astrocytes and pericytes of amygdala and hippocampus. A cell-level MAST model with
the animal as a random effect agrees: `DAG1` is detected in 81% of control astrocytes and 63% of infected ones
(80% → 59% in GFAP-high astrocytes). See [`macaque/README.md`](macaque/README.md) for the numbers, the caveats, and
why the conventional cell-level P-values are off by up to 291 orders of magnitude here.

## Repository layout

Each analysis is a self-contained folder with its own `README.md`, a rendered `index.html` report,
the scripts that produced it, and its figures. Further analyses will be added as sibling folders and
linked from the landing page.

```
.
├── index.html     landing page (GitHub Pages)
└── macaque/       GSE342497 re-analysis — report, scripts, figures
```

**Only code, figures and rendered reports are versioned.** No sequencing data, count matrices or
`.rds` objects: rebuild them by running each analysis's scripts against its public accession.

## GitHub Pages

The landing page and reports are served from `main` at the repository root. If Pages is not yet
enabled: **Settings → Pages → Build and deployment → Deploy from a branch → `main` / `/ (root)`**.
The `.nojekyll` file keeps the files served exactly as committed.

---

The preprint has not been peer reviewed. The analyses here are exploratory and
hypothesis-generating; sample sizes are small and the caveats in each folder's README are not
optional reading.
