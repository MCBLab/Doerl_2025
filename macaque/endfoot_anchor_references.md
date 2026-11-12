# Evidence that DAG1, DMD, UTRN, LARGE1, SNTA1 and DTNA build the astrocyte endfoot anchor

Primary literature establishing that each gene is structurally required for the dystrophin–glycoprotein
complex (DGC) that holds AQP4 and Kir4.1 at the perivascular astrocyte endfoot. Every entry was
verified against PubMed (PMID, journal, volume, pages and DOI checked).

## The architecture in one paragraph

In the perivascular endfoot membrane, **α-dystroglycan** (encoded by `DAG1`) binds laminin, agrin and
perlecan in the vascular basement membrane. Its binding depends on a sugar chain, *matriglycan*,
synthesised by **LARGE1**. The transmembrane **β-dystroglycan** subunit passes the link inward to
**dystrophin** (`DMD`, mainly the short glial isoform Dp71) or its paralogue **utrophin** (`UTRN`),
which scaffold **α-syntrophin** (`SNTA1`) and **α-dystrobrevin** (`DTNA`). α-Syntrophin's PDZ domain
binds the AQP4 C-terminus directly. Remove any link and AQP4 loses its perivascular polarisation while
total AQP4 protein can stay nearly unchanged — which is why a transcript-level fall in the anchor is
meaningful even when `AQP4` itself does not move.

```
basement membrane   LAMA/LAMB/LAMC · AGRN · HSPG2 · NID1 · COL4
                            │  (binds matriglycan, built by LARGE1)
                     α-dystroglycan  ─┐
                                      │  DAG1
                     β-dystroglycan  ─┘
                            │
                   Dp71 (DMD)  /  utrophin (UTRN)
                            │
              α-syntrophin (SNTA1) ── α-dystrobrevin (DTNA)
                            │
                    AQP4  ·  Kir4.1 (KCNJ10)
```

## By gene

### `DAG1` — dystroglycan (the core link)

| Reference | What it establishes |
|---|---|
| Guadagno E, Moukhles H. *Laminin-induced aggregation of the inwardly rectifying potassium channel, Kir4.1, and the water-permeable channel, AQP4, via a dystroglycan-containing complex in astrocytes.* **Glia** 2004;47:138–49. [PMID 15185393](https://pubmed.ncbi.nlm.nih.gov/15185393/) · [doi:10.1002/glia.20039](https://doi.org/10.1002/glia.20039) | Direct demonstration that laminin clusters Kir4.1 and AQP4 in astrocytes through a dystroglycan complex. |
| Noell S, et al. *Evidence for a role of dystroglycan regulating the membrane architecture of astroglial endfeet.* **Eur J Neurosci** 2011;33:2179–86. [PMID 21501259](https://pubmed.ncbi.nlm.nih.gov/21501259/) · [doi:10.1111/j.1460-9568.2011.07688.x](https://doi.org/10.1111/j.1460-9568.2011.07688.x) | Dystroglycan deletion changes endfoot membrane architecture and the orthogonal arrays of particles that AQP4 forms. |
| Satz JS, et al. *Distinct functions of glial and neuronal dystroglycan in the developing and adult mouse brain.* **J Neurosci** 2010;30:14560–72. [PMID 20980614](https://pubmed.ncbi.nlm.nih.gov/20980614/) · [doi:10.1523/JNEUROSCI.3247-10.2010](https://doi.org/10.1523/JNEUROSCI.3247-10.2010) | Glia-specific dystroglycan deletion separates the glial (basement-membrane/endfoot) role from the neuronal one. |
| Michele DE, et al. *Post-translational disruption of dystroglycan-ligand interactions in congenital muscular dystrophies.* **Nature** 2002;418:417–22. [PMID 12140558](https://pubmed.ncbi.nlm.nih.gov/12140558/) · [doi:10.1038/nature00837](https://doi.org/10.1038/nature00837) | Ligand binding requires the glycan, not just the protein — the result that makes `LARGE1` part of this module. |

### `DMD` — dystrophin, glial isoform Dp71

| Reference | What it establishes |
|---|---|
| Yang J, et al. *Dystrophin 71 deficiency causes impaired aquaporin-4 polarization contributing to glymphatic dysfunction and brain edema in cerebral ischemia.* **Neurobiol Dis** 2024;199:106586. [PMID 38950712](https://pubmed.ncbi.nlm.nih.gov/38950712/) · [doi:10.1016/j.nbd.2024.106586](https://doi.org/10.1016/j.nbd.2024.106586) | The most direct statement of the chain this analysis cares about: lose Dp71 → AQP4 depolarises → **glymphatic** function fails. |
| Frigeri A, et al. *Aquaporin-4 deficiency in skeletal muscle and brain of dystrophic mdx mice.* **FASEB J** 2001;15:90–8. [PMID 11149896](https://pubmed.ncbi.nlm.nih.gov/11149896/) · [doi:10.1096/fj.00-0260com](https://doi.org/10.1096/fj.00-0260com) | First demonstration that dystrophin loss reduces AQP4 at glial endfeet. |
| Vajda Z, et al. *Delayed onset of brain edema and mislocalization of aquaporin-4 in dystrophin-null transgenic mice.* **PNAS** 2002;99:13131–6. [PMID 12232046](https://pubmed.ncbi.nlm.nih.gov/12232046/) · [doi:10.1073/pnas.192457099](https://doi.org/10.1073/pnas.192457099) | Functional consequence: mislocalised AQP4 changes water handling, not just staining. |
| Nicchia GP, et al. *Dystrophin-dependent and -independent AQP4 pools are expressed in the mouse brain.* **Glia** 2008;56:869–76. [PMID 18338794](https://pubmed.ncbi.nlm.nih.gov/18338794/) · [doi:10.1002/glia.20661](https://doi.org/10.1002/glia.20661) | Only the perivascular pool is dystrophin-dependent — total AQP4 can look normal while polarisation is lost. |
| Nico B, et al. *Severe alterations of endothelial and glial cells in the blood-brain barrier of dystrophic mdx mice.* **Glia** 2003;42:235–51. [PMID 12673830](https://pubmed.ncbi.nlm.nih.gov/12673830/) · [doi:10.1002/glia.10216](https://doi.org/10.1002/glia.10216) | Extends the phenotype to barrier structure. |

### `UTRN` — utrophin (the paralogue that substitutes for dystrophin)

| Reference | What it establishes |
|---|---|
| Haenggi T, et al. *The role of utrophin and Dp71 for assembly of different dystrophin-associated protein complexes (DPCs) in the choroid plexus and microvasculature of the brain.* **Neuroscience** 2004;129:403–13. [PMID 15501597](https://pubmed.ncbi.nlm.nih.gov/15501597/) · [doi:10.1016/j.neuroscience.2004.06.079](https://doi.org/10.1016/j.neuroscience.2004.06.079) | Utrophin and Dp71 assemble *distinct* DGC variants at brain microvasculature — the reason both belong in the panel. |
| Haenggi T, Fritschy JM. *Role of dystrophin and utrophin for assembly and function of the dystrophin glycoprotein complex in non-muscle tissue.* **Cell Mol Life Sci** 2006;63:1614–31. [PMID 16710609](https://pubmed.ncbi.nlm.nih.gov/16710609/) · [doi:10.1007/s00018-005-5461-0](https://doi.org/10.1007/s00018-005-5461-0) | The standard review of the brain DGC; good single citation for the architecture as a whole. |
| Blake DJ, et al. *Different dystrophin-like complexes are expressed in neurons and glia.* **J Cell Biol** 1999;147:645–58. [PMID 10545507](https://pubmed.ncbi.nlm.nih.gov/10545507/) · [doi:10.1083/jcb.147.3.645](https://doi.org/10.1083/jcb.147.3.645) | Establishes that glia run their own DGC composition, distinct from neurons. |

### `LARGE1` — the glycosyltransferase that makes the link functional

| Reference | What it establishes |
|---|---|
| Michele DE, et al. **Nature** 2002;418:417–22. [PMID 12140558](https://pubmed.ncbi.nlm.nih.gov/12140558/) · [doi:10.1038/nature00837](https://doi.org/10.1038/nature00837) | Hypoglycosylated α-dystroglycan cannot bind laminin or agrin: the anchor fails without the sugar. |
| Inamori K, et al. *Dystroglycan function requires xylosyl- and glucuronyltransferase activities of LARGE.* **Science** 2012;335:93–6. [PMID 22223806](https://pubmed.ncbi.nlm.nih.gov/22223806/) · [doi:10.1126/science.1214115](https://doi.org/10.1126/science.1214115) | Identifies LARGE as the enzyme building matriglycan, the actual ligand-binding moiety. |
| Zhou M, et al. *Large is required for normal astrocyte migration and retinal vasculature development.* **Cell Biosci** 2017;7:18. [PMID 28428837](https://pubmed.ncbi.nlm.nih.gov/28428837/) · [doi:10.1186/s13578-017-0143-9](https://doi.org/10.1186/s13578-017-0143-9) | LARGE loss disturbs astrocyte–vasculature interaction specifically. |

### `SNTA1` — α-syntrophin (binds AQP4 directly)

| Reference | What it establishes |
|---|---|
| Neely JD, et al. *Syntrophin-dependent expression and localization of Aquaporin-4 water channel protein.* **PNAS** 2001;98:14108–13. [PMID 11717465](https://pubmed.ncbi.nlm.nih.gov/11717465/) · [doi:10.1073/pnas.241508198](https://doi.org/10.1073/pnas.241508198) | The founding result: α-syntrophin is required for AQP4 to reach and stay at the endfoot. |
| Amiry-Moghaddam M, et al. *An α-syntrophin-dependent pool of AQP4 in astroglial end-feet confers bidirectional water flow between blood and brain.* **PNAS** 2003;100:2106–11. [PMID 12578959](https://pubmed.ncbi.nlm.nih.gov/12578959/) · [doi:10.1073/pnas.0437946100](https://doi.org/10.1073/pnas.0437946100) | Defines the perivascular AQP4 pool functionally. |
| Amiry-Moghaddam M, et al. *Delayed K⁺ clearance associated with aquaporin-4 mislocalization: phenotypic defects in brains of α-syntrophin-null mice.* **PNAS** 2003;100:13615–20. [PMID 14597704](https://pubmed.ncbi.nlm.nih.gov/14597704/) · [doi:10.1073/pnas.2336064100](https://doi.org/10.1073/pnas.2336064100) | Couples the anchor to K⁺ homeostasis — why `KCNJ10` moves with it. |
| Bragg AD, et al. *Assembly of a perivascular astrocyte protein scaffold at the mammalian blood-brain barrier is dependent on α-syntrophin.* **Glia** 2006;53:879–90. [PMID 16609960](https://pubmed.ncbi.nlm.nih.gov/16609960/) · [doi:10.1002/glia.20347](https://doi.org/10.1002/glia.20347) | Shows the whole perivascular scaffold — syntrophin, dystrobrevin, AQP4 — assembles as one unit. |

### `DTNA` — α-dystrobrevin

| Reference | What it establishes |
|---|---|
| Lien CF, et al. *Absence of glial α-dystrobrevin causes abnormalities of the blood-brain barrier and progressive brain edema.* **J Biol Chem** 2012;287:41374–85. [PMID 23043099](https://pubmed.ncbi.nlm.nih.gov/23043099/) · [doi:10.1074/jbc.M112.400044](https://doi.org/10.1074/jbc.M112.400044) | Glia-specific α-dystrobrevin loss → AQP4 redistribution, barrier failure, progressive oedema. |
| Bragg AD, et al. **Glia** 2006;53:879–90. [PMID 16609960](https://pubmed.ncbi.nlm.nih.gov/16609960/) · [doi:10.1002/glia.20347](https://doi.org/10.1002/glia.20347) | Places α-dystrobrevin in the same perivascular scaffold as α-syntrophin and AQP4. |
| Blake DJ, et al. **J Cell Biol** 1999;147:645–58. [PMID 10545507](https://pubmed.ncbi.nlm.nih.gov/10545507/) · [doi:10.1083/jcb.147.3.645](https://doi.org/10.1083/jcb.147.3.645) | Characterises the dystrobrevin-containing glial complex. |

## Why this matters for the ZIKV comparison

[Doerl et al. 2025](https://doi.org/10.1101/2025.10.26.684625) report reduced AQP4 coverage of vessels
and abnormal AQP4 patterns when ZIKV infects before glymphatic maturation. The literature above says
that exactly this phenotype — AQP4 present but depolarised — is what happens when the anchor fails,
and Yang et al. 2024 closes the loop by showing Dp71 loss impairs glymphatic function specifically. In
the macaque data, `AQP4` transcript is unchanged while `DAG1`, `LARGE1`, `UTRN` and the
basement-membrane ligands fall, which is the transcriptional signature this architecture predicts for
a loss of polarisation rather than a loss of channel.

**This remains an inference.** None of it is measured here: mRNA is not protein and certainly not
localisation, and the macaque cohort is 3 controls versus 2 infected animals. AQP4 and
β-dystroglycan immunostaining scored as perivascular polarisation is the experiment that would test it.
