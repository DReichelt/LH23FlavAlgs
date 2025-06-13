# Comparing flavoured jet algorithms | PhysTeV 2023

Analyses as used in [arXiv:2506.XXXXX](https://arxiv.org/abs/2506.XXXXX) to compare apporaches to infrared save definition of anti - $k_t$ jets.

References for algorithms:

- SDF: [arXiv:2205.01109](https://arxiv.org/abs/2205.01109)
- CMP: [arXiv:2205.11879](https://arxiv.org/abs/2205.11879) 
- GHS: [arXiv:2208.11138](https://arxiv.org/abs/2208.11138)
- IFN: [arXiv:2306.07314](https://arxiv.org/abs/2306.07314)
- WTA: [arXiv:2205.01117](https://arxiv.org/abs/2205.01117)

Code implementing the algorithms is available as [fastjet](https://fastjet.fr) contribs.

## Details on analyses

The analyses can be compiled as plugins to the [Rivet](https://rivet.hepforge.org) analysis framework, see documentation there for details, and depend on having fastjet contrib installed.

### FlavAlgAnalysis

Analysis used for plots in Appendix C of [arXiv:2506.XXXXX](https://arxiv.org/abs/2506.XXXXX). 

The selection criteria are loosely inspired by the CMS analysis [arXiv:1611.06507](https://arxiv.org/abs/1611.06507), and the corresponding [Rivet analysis code](https://rivet.hepforge.org/analyses/CMS_2017_I1499471.html). `FlavAlgAnalysis.yoda` contains the original reference files from Rivet, used to inherit binning for some histograms.
