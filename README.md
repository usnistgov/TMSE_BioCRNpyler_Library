# Toehold-mediated strand exchange (TMSE) circuits 

![TMSE-BioCRNpyler Library](../schematics/Figure_readme.png)

## Description

Toehold-mediated strand exchange circuits (TMSE)-BioCRNpyler library is an open source package library to [BioCRNpyler][BioCRNPyler]. 
The library can be used to model toehold-mediated strand exchange circuits as we have demonstrated in [TMSE-paper][TMSE_doi]. 

The new implemented components and mechanisms included in this library are listed below:

**Components**
- `RNAassembly`: High-level representation of an RNA expression construct.
- `LeakctRSD`: Gate DNA that exhibit leak behavior, allowing for unintended transcription.
- `RNAmaturation`: RNA molecules that require processing prior to strand-exchange.
- `StrandExchangeMol`: RNA molecules that can undergo strand exchange.
- `TMSE`: Molecules that undergo two-domain mediated strand exchange.

**Mechanism**
- `MaturationO1`: First order maturation mechanism for preprocessed molecules.
- `GateTxLeak`: Leakage mechanism for unintended transcription of Gate DNA to Output strands.
- TMSE: Mechanism for strand exchange.
   - `OneStepTMSE`: One step catalytic mechanism for strand exchange.
   - `TwoStepTMSE`: Two step catalytic mechanism for strand exchange.
   - `ThreeStepTMSE`: Three step catalytic mechanism for strand exchange.

----

## Getting started / installation

To get started using BioCRNpyler to model TMSE circuits we suggested installing BioCRNpyler following the [Installation Guide][install-guide]. Once BioCRNpyler is installed to get started using the TMSE package by:
```
git clone https://github.com/usnistgov/TMSE_BioCRNpyler_Library
```

### Dependencies

- [BioCRNpyler](https://github.com/BuildACell/bioCRNpyler)
- [Bioscrape](https://github.com/biocircuits/bioscrape) 

## Related Material

Associated project pages on the NIST website:
   - [RNA Computation and Metrology for Engineering Biology][ctRSD-project]
   <!-- - [Website 2][TMSEmodel-project] -->

## How to cite

An article about the application of the extended classes is available in bioRxiv. 
If this BioCRNpyler Library helped you in your research, please cite:

```
@article{tmse_biocrnpyler2026,
  title={A modeling pipeline for toehold-mediated strand exchange circuits across synthetic biology applications},
  author={Jurado, Zoila and Taghon, Geoffrey and Schaffter, Samuel W.},
  journal = {bioRxiv},
  year = {2026},
  doi = {10.64898/2026.07.22.739222},
  url = {https://www.biorxiv.org/content/10.64898/2026.07.22.739222v1},
  note = {Available at: \url{https://github.com/usnistgov/TMSE_BioCRNpyler_Library}}
}
```

or the GitHub site: https://github.com/usnistgov/TMSE_BioCRNpyler_Library. 

An [article][BioCRNPyler_paper] about BioCRNpyler is available in PLOS Computational Biology. Please also site the orignal BioCRNpyler package:

```
@article{biocrnpyler2022,
  title={{BioCRNpyler}: {Compiling} chemical reaction networks from
         biomolecular parts in diverse contexts},
  author={W. Poole and A. Pandey and A. Shur and Z. A. Tuza and
          R. M. Murray},
  journal = {PLOS Computational Biology},
  year = {2022},
  month = {04},
  volume = {18},
  url = {https://doi.org/10.1371/journal.pcbi.1009987},
  pages = {1-19}
}
```

or the GitHub site: https://github.com/BuildACell/BioCRNPyler.

## Contact

- [Zoila Jurado][zj-nist]
- [Sam Schaffter][sws-nist]

## License

See [LICENSE.md](LICENSE.md). This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to Title 17 Section 105 of the United States Code, this software is not subject to copyright protection and is in the public domain.

<!-- References -->
[BioCRNPyler]: https://github.com/BuildACell/BioCRNPyler
[install-guide]: https://github.com/BuildACell/BioCRNPyler?tab=readme-ov-file#installation
[BioCRNPyler_paper]: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009987
[ctRSD-project]: https://www.nist.gov/programs-projects/rna-computation-and-metrology-engineering-biology
[cellenggroup]: [https://www.nist.gov/mml/bbd/cellular-engineering-group]
[sws-nist]:{https://www.nist.gov/people/samuel-schaffter}
[zj-nist]:{https://www.nist.gov/people/zoila-jurado}
[TMSE_doi]: https://www.biorxiv.org/content/10.64898/2026.07.22.739222v1
<!-- [TMSEmodel-project]: -->
