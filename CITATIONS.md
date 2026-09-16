# nf-core/crispr: Citations

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

### Statistical methods behind the inference step

The pipeline's two inference methods both test conditional independence by resampling a
perturbation's assignment. The conditional randomization test (CRT) they rest on, its power
properties, its application to single-cell CRISPR screens, and the saddlepoint approximation
that makes it affordable at screen scale are due to:

> Candès E, Fan Y, Janson L, Lv J. Panning for gold: 'model-X' knockoffs for high dimensional controlled variable selection. J R Stat Soc Series B Stat Methodol. 2018;80(3):551-577. doi: 10.1111/rssb.12265. (Introduces the conditional randomization test.)

> Katsevich E, Ramdas A. On the power of conditional independence testing under model-X. Electron J Stat. 2022;16(2):6348-6394. doi: 10.1214/22-EJS2085.

> Barry T, Wang X, Morris JA, Roeder K, Katsevich E. SCEPTRE improves calibration and sensitivity in single-cell CRISPR screen analysis. Genome Biol. 2021;22. doi: 10.1186/s13059-021-02545-2. (SCEPTRE; the CRT applied to single-cell CRISPR screens.)

> Barry T, Mason K, Roeder K, Katsevich E. Robust differential expression testing for single-cell CRISPR screens at low multiplicity of infection. Genome Biol. 2024;25. doi: 10.1186/s13059-024-03254-2. (Low-MOI analysis, the effective-sample-size diagnostic, and the non-targeting-cell contrast this pipeline's `INFERENCE_control_group = 'nt_cells'` selects.)

> Niu Z, Huang Z, Ray Choudhury J, Katsevich E. Saddlepoint approximations for plug-in resampling. arXiv:2407.08911. (Introduces spaCRT, the saddlepoint approximation to the distilled CRT.)

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
