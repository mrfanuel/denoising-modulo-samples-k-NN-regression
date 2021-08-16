# denoising-modulo-samples-k-NN-regression

This repository contains scripts associated to the paper:
>M. Fanuel and H. Tyagi, Denoising modulo samples: k-NN regression and tightness of SDP relaxation, accepted in Information and Inference: A Journal of the IMA
https://arxiv.org/pdf/2009.04850.pdf

A few implementations of denoising/unwrapping methods were shared with us by M. Cucuringu; see the following paper
>M. Cucuringu and H. Tyagi. Provably robust estimation of modulo 1 samples of a smooth function with applications to phase unwrapping. Journal of Machine Learning Research, 21(32):1–77, 2020.
https://www.jmlr.org/papers/volume21/18-143/18-143.pdf

## files

- Denoise_1D_main_ex1.m reproduces Figure 2 and Figure 3.
- Denoise_1D_main_ex2.m reproduces Figure 4 and Figure 5.
- Denoise_2D_Vesuvius.m reproduces Figure 6.

### Vesuvius data

The elevation map of MountVesuvius (N40E014.hgt.zip) was downloaded thanks to the readhgt.m script written by [François Beauducel](https://www.ipgp.fr/~beaudu/); see https://dds.cr.usgs.gov/srtm/version2_1.
This data was stored in Vesuvius_zoom.mat file after a preprocessing.

## Dependencies

These scripts require Matlab Statistics and Machine Learning Toolbox (version Matlab R2019b).



