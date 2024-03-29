---
title: 'SeisMix: A Python package for estimating mixing rates from seismic oceanographic images'
tags:
  - Python
  - oceanography
  - turbulence
  - internal waves
  - seismic
authors:
  - name: Alex Dickinson
    orcid: 0000-0001-7184-927X
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Bullard Laboratories, Department of Earth Sciences, University of Cambridge, Cambridge, United Kingdom
   index: 1
 - name: Now at School of Natural and Environmental Sciences, Newcastle University, Newcastle-upon-Tyne, United Kingdom
   index: 2
date: \? \? 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Marine seismic reflection data provide an unprecedented way to image changes in temperature and salinity within the oceanic water column.  

@Dickinson2017

@Dickinson2020

@Gunn2021

# Statement of Need
As of ? 2021, (25 in Dickinson and Gunn + more) English-language articles have presented results based on spectral analysis of seismic oceanographic images / tracked seismic reflections (see Dickinson2022 for a review). The majority of these articles  Each group of researchers has developed its own software, and so it has been difficult to compare results from different studies. Many analyses have relied in part on closed-source, proprietary software that is not available to all researchers. 


<!-- However, there exists no standard approach to spectral analysis, and so it has been difficult to compare results from different studies. Many analyses have relied on closed-source, proprietary software that is not available to all researchers.  -->



`SeisMix` is an open-source Python package that provides a fully automated workflow for spectral analysis of seismic oceanographic images. It assigns rigorous uncertainties to all estimated values. Using the standardised approach of `SeisMix`, results from different datasets can be easily compared.


 <!-- can standardise analysis of different datasets. -->

The scripts underlying `SeisMix` have produced results for three published research papers (@Dickinson2017, @Dickinson2020, @Gunn2021).

# Documentation and examples

`SeisMix` contains modules for:
- Reading input data in [SEG-Y](https://wiki.seg.org/wiki/SEG-Y) format
- Assessing the spectral content of a seismic image
- Tracking seismic reflections across an image
- Computing horizontal-wavenumber spectra from tracked reflections
- Estimating diapycnal diffusivity from computed spectra


`SeisMix` is provided with [Jupyter](https://jupyter.org/) notebooks that demonstrate its use. The notebooks cover:

- Estimation of diapycnal diffusivity from the spectral signal of turbulence

- Estimation of diapyncal diffusivity from the spectral signal of internal waves

- Fitting of the Garrett-Munk model spectrum to observed internal-wave spectra


I hope that `SeisMix` will be of use to the seismic oceanographic community and will encourage development of other standardised methods of analysis.
Please report any errors to nad38@cantab.ac.uk.


<!-- Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% } -->

# Acknowledgements

Thanks to [Kathy Gunn](https://kathygunn.github.io/) for testing `SeisMix` and to [Bryn Pickering](https://www.brynpickering.com/) for helping build the distribution. Thanks also to [Jody Klymak](https://ocean-physics.seos.uvic.ca/~jklymak/) for generously providing access to his [MATLAB toolbox](jklymak.github.io/GarrettMunkMatlab) for calculating Garrett-Munk spectra. Parts of `SeisMix` were inspired by the work of [Katy Sheen](https://geography.exeter.ac.uk/staff/index.php?web_id=Katy_Sheen) and Matthew Falder. Development of `SeisMix` was funded by the [Natural Environment Research Council](https://nerc.ukri.org/) and by the [North East Local Enterprise Partnership](https://www.northeastlep.co.uk/).

# References