<!---[![Chat on Gitter](https://img.shields.io/gitter/room/calliope-project/calliope.svg?style=flat-square)](https://gitter.im/calliope-project/calliope)--->

<!-- [![Master branch build status](https://img.shields.io/azure-devops/build/calliope-project/371cbbaa-fa6b-4efb-9b23-c4283a8e33eb/1?style=flat-square)](https://dev.azure.com/calliope-project/calliope/_build?definitionId=1) -->
<!-- [![Documentation build status](https://img.shields.io/readthedocs/calliope.svg?style=flat-square)](https://readthedocs.org/projects/calliope/builds/)
[![Test coverage](https://img.shields.io/codecov/c/github/calliope-project/calliope?style=flat-square&token=b4fd170f0e7b43679a8bf649719e1cea)](https://codecov.io/gh/calliope-project/calliope)
[![PyPI version](https://img.shields.io/pypi/v/calliope.svg?style=flat-square)](https://pypi.python.org/pypi/calliope)
[![Anaconda.org/conda-forge version](https://img.shields.io/conda/vn/conda-forge/calliope.svg?style=flat-square&label=conda)](https://anaconda.org/conda-forge/calliope)
[![JOSS DOI](https://img.shields.io/badge/JOSS-10.21105/joss.00825-green.svg?style=flat-square)](https://doi.org/10.21105/joss.00825) -->


<!-- # Mixing codes for Kathy -->

<!-- Make sure code works with example data from FSC -->

<!-- Unused modules?:
contours_to_midpoints_fs3d.py
contours_to_midpoints.py -->


Base dependencies:

Python 3.8.12

- [`matplotlib >= 3.3.3`](https://matplotlib.org/)
- [`numpy >= 1.19.0`](http://numpy.org)
- [`openpyxl >= 3.0.9`](https://openpyxl.readthedocs.io/en/stable/)
- [`pandas >= 1.0.5`](https://pandas.pydata.org/)
- [`scipy >= 1.5.1`](https://scipy.org)


Extra dependencies for running example notebooks:

- [`jupyter`](https://jupyter-notebook.readthedocs.io/en/stable/)
- [`jupyterlab`](https://jupyterlab.readthedocs.io/en/stable/) 

for running notebooks using either Jupyter Notebook or JupyterLab.


To set up SeisMix, first create conda environment:

```conda create --name NAME```

Activate conda environment:

```conda activate NAME```

Install Python 3.8 and jupyter:

```conda install python=3.8 jupyter```



To install SeisMix with its base dependencies:

Navigate to directory SeisMix



```bash
pip install .
```

To install extra dependencies for running example notebooks:

```bash
pip install .[examples]
```


Installation with pip:

Create a virtual environment (called env_SeisMix in this example):


```bash 
pip install .
```

To install optional packages for running example jupyter notebooks:
```bash 
pip install .[examples]
```

To create a Jupyter kernel for running the notebooks, use
ipython kernel install --name "local-venv" --user


Needs packages:
segyio: https://pypi.org/project/segyio/ - this is currently used only in example notebook for reading data in. Build function for reading segy in

Dependencies:
numpy
numpy.matlib

scipy.fftpack
scipy.signal
scipy.stats
scipy.optimize
scipy.interpolate

mtspec

matplotlib.pyplot

gsw


To add: functionality to read CTD profile for density etc

``numpy``

`test`

## Installation

<!-- ### Using pip

Set up virtual pip environment.

```bash
pip install SeisMix -r requirements.txt
``` -->



### Dependencies

SeisMix uses the following packages:

- [`mtspec`](https://krischer.github.io/mtspec/)
- [`numpy`](http://numpy.org)
- [`scipy`](https://scipy.org)
- [`segyio`](https://segyio.readthedocs.io/en/1.5.3/index.html)

__Optional dependencies__ for running the example Notebooks:
- [`gsw`](https://github.com/TEOS-10/GSW-python)
- [`jupyter`](https://jupyter.org/)
- [`matplotlib`](https://matplotlib.org/)



# Maths

All text currently copied from papers - need to vary

## Noise analysis

### Estimation of signal-to-noise ratio

`noise_analysis.estimate_signal_to_noise`

analyze the frequency content of the seismic image  by defining the signal-to-noise ratio between two adjacent traces as

$$ \textrm{SNR} = \sqrt{\frac{|c|}{|a - c|}}, $$

where $c$ is the maximum value of the cross-correlation between these traces,  and $a$ is the value of the zero-lag autocorrelation of the first trace @Holbrook2013. The signal-to-noise ratio for a given panel  is gauged using  the median value of signal-to-noise ratios for  adjacent traces.

### Computation of amplitude and phase spectra

`noise_analysis.make_amplitude_and_phase_spectra`

### Computation of direct data transform

`noise_analysis.make_direct_data_transform`

## Multitaper Fourier transform

Multi-taper Fourier transforms are computed using the `mtspec` package.
This package is used in the functions:
- `noise_analysis.make_amplitude_and_phase_spectra`
- 

Estimation of multi-taper Fourier transforms is controlled by two parameters: $j_{res}$ (related to the wavenumber resolution: increasing $j_{res}$ decreases wavenumber resolution) and $K_{res}$ (related to the power resolution:  increasing $K_{res}$ increases power resolution). $K_{res} = 2j_{res}-2$ represents the optimal power resolution for a given value of $j_{res}$. We choose 
$j_{res}=4, K_{res}=6$ as an appropriate  compromise between wavenumber resolution and power resolution. The wavenumber resolution can be approximately estimated using

$$ j_{res} \frac{2k_N}{N} = \frac{j_{res}}{N \Delta x}, $$

where $k_N$ is the Nyquist wavenumber, $N$ is the number of data points in the input signal, and $\Delta x$ is the cmp spacing.

Following @Percival1993, the width of the $100(1-2p) \%$ confidence interval for $\log_{10}\left(\phi_{\tau}(k_x)\right)$ is given by

$$ \log_{10}(e) \log_{10}\left( \frac{Q_{2K_{res}}(1-p)}{Q_{2K_{res}}(p)} \right), $$

where $Q_{2K_{res}}(p)$ is the $p^{th}$ percentage point for the $\chi^2$ distribution with $2K_{res}$ degrees of freedom and $K_{res}$ is the multi-taper parameter described above.



## Depth conversion

$\phi_{\tau}(k_x)$ is converted to a displacement power spectrum expressed in terms of depth, $\phi_{\xi}(k_x)$, by

$$ \phi_{\xi}(k_x) = \frac{v^2}{4} \phi_{\tau}(k_x), $$

where $v$ is the value of a regionally averaged profile of sound speed at the average depth of the tracked reflection. The error in $\log_{10}\left(\phi_\xi(k_x)\right)$, $\sigma_{\phi_\xi}$, is given by

$$ \sigma_{\phi_\xi}^2 = \sigma_{\phi_\tau}^2 + \left(2 \log_{10}(e) \frac{\sigma_v}{v} \right)^2,$$

where $\sigma_v$ is the standard deviation of $v$. An error of $\pm$30~m~s$^{-1}$ in sound speed corresponds to $\frac{\sigma_v}{v} \sim 0.02$, and so errors in sound speed have an insignificant effect.

### Conversion to horizontal-gradient spectra

Displacement power spectra, $\phi_{\xi}(k_x)$, are converted to spectra of the horizontal gradient of vertical displacement, $\phi_{\xi_x}(k_x)$, using 

$$ \phi_{\xi_x}(k_x) = \left( 2 \pi k_x \right)^2 \phi_{\xi}(k_x) \log_{10}\left(\phi_{\xi_x}(k_x)\right) = 2 \log_{10} \left( 2\pi k_x \right) + \log_{10}\left(\phi_{\xi}(k_x)\right) $$

Ignoring the effects of spectral leakage, the error in $\log_{10}\left(\phi_{\xi_x}(k_x)\right)$ is  equal to $\sigma_{\phi_\xi}$.

## Spectral fitting
Rewrite text as from Dickinson (2020): Four parameters control the form of the model: the gradient, $a$, of the low-wavenumber subrange; the coordinates, $(b,c)$, of the intersection point between the low-wavenumber and turbulent regimes; and the width, $d$, of the turbulent subrange. The model is fitted by defining a misfit function 

$$ M_1 = \sum_{i=0}^M \left\[\left\[\log_{10}(\phi_{\xi_x}^m) \right\]_{i} - \left\[\log_{10}( \left\langle \phi_{\xi_x} \right\rangle_p ) \right]_{i} \right\]^2 $$

$$M\_1 = \\sum\_{i=1}^D \\left\[\\left\[\\log\_{10}(\\phi\_{\\xi_x}^m) \\right\]\_i - \\left\[\\log\_{10}( \\left\\langle \\phi\_{\\xi_x} \\right\\rangle_p ) \\right\]\_i\\right\]^2,$$


where $\phi_{\xi_x}^m$ is the model and $M$ is the number of points in the power spectrum. Best-fitting values of $a$, $b$, $c$ and $d$ are found by minimizing $M_1$ using Powell's conjugate direction algorithm @Powell1964.






## Internal-wave subrange
### Fitting of straight lines

Measurements by towed instruments suggest that, for $10^{-3} \lesssim k_x \lesssim 10^{-2}$ cpm, horizontal wavenumber spectra are reasonably approximated by a power law that takes the  form

$$ \phi_{\xi_x} \propto k_x^p, $$

where $p$ is within the range of  $-0.5$ to $0$ \citep[e.g.,][]{Katz1973, McKean1974, Zenk1975, Bell1976, Muller1978, Katz1979}. 
Based on such   observations, a revised version  of the {\sf GM72} model proposed by \citet[][henceforth {\sf GM75}]{Garrett1975} exploits   a power law description where  $p = -0.5$ at  the limit of high horizontal wavenumber (i.e., $k_x \gtrsim 10^{-2}$ cpm).

These low-wavenumber subranges are fitted with straight lines using least squares linear regression. A  standard error, $\sigma_{\hat{m}}$, on the fitted gradient, $\hat{m}$, is estimated from the covariance matrix. $\sigma_{\phi_\xi}$ is taken as the error on the individual data points. The weighted mean, $\bar{m}_w $, of the fitted gradients is given by

$$ \bar{m}\_w = \frac{\sum\limits_{i=0}^M \hat{m}_i (\sigma_{\phi_\xi})\_i^{-2}}{\sum\limits_{i=0}^M (\sigma_{\phi_\xi})\_i^{-2}}, $$

where $M$ is the total number of low-wavenumber subranges.
Wrong: We note that there is no accepted definition of the standard error on the weighted mean, and so we quote the standard deviation of the sample as the error.




## Estimation of diapycnal diffusivity

### Estimation of turbulent dissipation rate and conversion to diapycnal diffusivity

The diapycnal diffusivity, $K$, is estimated from $\epsilon$ using the Osborn relationship  so that
\label{eq:gom_osborn_relation}

$$ K =  \frac{\frac{g}{\rho}\langle w'\rho' \rangle}{\frac{-g}{\rho}\frac{d\rho}{dz}} = \frac{B}{\bar{N}^2} = \frac{\Gamma \epsilon}{\bar{N}^2}, $$

where $g$ is gravitational  acceleration, $\rho$ is  mean density, $\rho'$ is  density perturbation, $w'$ is  vertical component of the velocity perturbation, $B$ is  buoyancy flux, and $\Gamma$ is the dissipation flux coefficient \citep{Osborn1980}. Angular brackets denote averaging over time. This  method assumes a balance between turbulent production, dissipation, and work done against buoyancy for the case of a statistically steady flow. Following convention, we assume a constant dissipation flux coefficient of $\Gamma = 0.2$, which represents  the upper bound proposed by \citet{Osborn1980}. It is now recognized  that the assumption of  constant $\Gamma$ is undoubtedly a significant  simplification \citep[e.g.,][]{Mashayek2017}.

### Using internal-wave subrange

Estimates of $\epsilon$ are determined  by reference to the dissipation rate, $\epsilon_r$, of a GM spectrum at a specified buoyancy frequency, $N_r$, and to a given  latitude with a Coriolis parameter $f_r$. Here, we use
\label{eq:gom_polzin_epsilon_parametrization}

$$ \epsilon = \epsilon_r \left( \frac{\bar{N}}{N_r} \right)^2 \frac{f \ \mathrm{arccosh} \left( \bar{N} / f \right)}{f_r \ \mathrm{arccosh} \left( N_r / f_r \right)} \frac{R_{\omega_{GM}}}{R_{\omega}} \frac{(R_\omega + 1)}{(R_{\omega_{GM}} + 1)} \sqrt{\frac{R_{\omega_{GM}} - 1}{R_\omega - 1}} \hat{E}^2, $$

where $\bar{N}$ is the local mean stable stratification, $f$ is the local Coriolis parameter, $R_\omega$ is the local shear-strain ratio, $R_{\omega_{GM}} = 3$ cph is the shear-strain ratio of the GM wave field, and $\hat{E}$ is a measure of the energy of the observed wave field \citep{Polzin2014}. $\hat{E}$ is evaluated in terms of horizontal gradient spectra using 
\label{eq:energy_level_slope}

$$ \hat{E} = \frac{R_\omega}{R_{\omega_{GM}}} \hat{\phi}\_{\xi_x} = \frac{R_\omega}{R_{\omega_{GM}}} \frac{\left\langle \phi_{\xi_x} \right\rangle}{\left\langle \phi_{\xi_x}^{GM} \right\rangle}, $$

where $\phi_{\xi_x}^{GM}$ is the  GM spectrum for the  local value of the buoyancy frequency, $\bar{N}$. Angular brackets denote  integration over the wavenumber range of the internal-wave subrange.

Following \citet{Gregg2003}, we compare observed dissipation rates with  those predicted by a GM model at 30\textdegree \ latitude. We choose {\sf GM76} as the reference spectrum with a mode number scale $j_\star = 4$ and a high wavenumber saturated subrange in the vertical spectra for $k_z > 0.1$ cpm. Equation \ref{eq:gom_polzin_epsilon_parametrization} now  becomes
 \label{eq:gom_polzin_epsilon_gm76_parametrization}
 
 $$ \epsilon =  7.2 \times 10^{-10}  \ \left( \frac{\bar{N}}{N_r} \right)^2 \frac{f \ \mathrm{arccosh} \left( \bar{N} / f \right)}{f_r \ \mathrm{arccosh} \left( N_r / f_r \right)} \frac{R_{\omega} (R_\omega + 1)}{12} \sqrt{\frac{2}{R_\omega - 1}} \frac{\left\langle \phi_{\xi_x} \right\rangle^2}{\left\langle \phi_{\xi_{x}}^{GM} \right\rangle^2}, $$

with  units of $\mathrm{m}^2$ $\mathrm{s}^{-3}$.

$\epsilon$ and $K$ are estimated from internal wave spectral subranges using Equations \ref{eq:gom_polzin_epsilon_parametrization}, \ref{eq:energy_level_slope} and \ref{eq:gom_osborn_relation}. Combination of these equations shows that

$$ K = K_r L(f, \bar{N}) J(R_\omega) \frac{\left\langle \phi_{\xi_x} \right\rangle^2}{\left\langle \phi_{\xi_x}^{GM} \right\rangle^2}, $$
\label{eq:gom_polzin_k_parametrization}

where

$$ K_r = \frac{\Gamma \epsilon_r}{N_r^2}, $$

$$ L(f, \bar{N}) = \frac{f \ \mathrm{arccosh} \left( \bar{N} / f \right)}{f_r \ \mathrm{arccosh} \left( N_r / f_r \right)} $$

and

$$ J(R_\omega) = \frac{1}{6 \sqrt{2}} \frac{R_\omega(R_\omega + 1)}{\sqrt{R_\omega - 1}}. $$

Estimated values of $K$ have only weak explicit dependences on $\bar{N}$ and $f$ through the term $L(f, \bar{N})$.
Estimated values of $K$ also have an implicit dependence on $\bar{N}$ through $\phi_{\xi_x}^{GM}$, which is approximately proportional to $\bar{N}^{-1}$. Uncertainty in $\bar{N}$ translates to an approximate error in $\log_{10}(K)$ of

$$ \left| \frac{\partial}{\partial \bar{N}} \left[ \log_{10}\left(L(f, \bar{N})\right) + 2 \log_{10}(\bar{N}) \right] \right| \sigma(\bar{N}) \\
 = \log_{10}(e) \left| \frac{1}{\ \mathrm{arccosh}(\bar{N}/f)}  \frac{1}{f}  \frac{1}{\sqrt{(\bar{N}/f)^2 - 1}} + \frac{2}{\bar{N}} \right|  \sigma(\bar{N}), $$
 
where vertical bars indicate that the modulus is taken and  $\sigma(\bar{N})$ is the error in $\bar{N}$.

When gauging  the effect of uncertainties in $R_\omega$ and $\Gamma$, we have not formally propagated errors since  the underlying distributions of $R_\omega$ and $\Gamma$ are not known with sufficient accuracy. Instead, we consider the effect of choosing upper and lower bounds for $R_\omega$ and $\Gamma$ on $K$.

We have no means to estimate $R_\omega$, which has  an average value of $7 \pm 3$ within the abyssal ocean \citep{Kunze2006}. It is plausible that $R_\omega$  varies systematically when   water depth  shallows toward  the continental slope. Even if $R_\omega$ varies between 2 and 20, the  assumption of a constant value of $R_\omega = 7$ yields  an uncertainty  in $\log_{10}(K)$ of 
$<0.6$ logarithmic units such that

$$ \log_{10}(K) = \log_{10}\left(6\sqrt{2} J(R_\omega)\right) +  A \\
 = \log_{10} \left( \frac{R_\omega(R_\omega + 1)}{\sqrt{R_\omega - 1}} \right) + A \\
 \approx 0.8 + A \ \mathrm{for} \ R_\omega = 2 \\
 \approx 1.4 + A \ \mathrm{for} \ R_\omega = 7 \\
 \approx 2.0 + A \ \mathrm{for} \ R_\omega = 20, $$

where $A$ includes all terms independent of $R_\omega$. To account for all the observed spatial variation in $\log_{10}(K)$, $R_\omega$ would have to vary over three orders of magnitude.

Finally, $\Gamma$ is probably  not  constant but is likely to vary both temporally and spatially \citep{Mashayek2017}. Taking lower and upper bounds of $\Gamma \sim 0.1$ and $\Gamma \sim 0.4$, respectively, implies that assumption of a constant value of $\Gamma = 0.2$ introduces an uncertainty  in $\log_{10}(K)$ of up to 0.3 logarithmic units such that 

$$ \log_{10}(K) = \log_{10}(\Gamma) +  B \\
 = -1 + B \ \mathrm{for} \ \Gamma = 0.1 \\
 \approx -0.7 + B \ \mathrm{for} \ \Gamma = 0.2 \\
 \approx -0.4 + B \ \mathrm{for} \ \Gamma = 0.4, $$

where $B$ includes all terms independent of $\Gamma$. We also note that a constant flux Richardson number, $R_f$, of 0.17 is assumed in estimation of $\epsilon_r = 7.2 \times 10^{-10}$~m$^2$~s$^{-1}$. 

Error in the energy level ratio of individual spectra is estimated as the standard deviation of the ratio

$$ \frac{\phi_{\xi_x}}{\phi_{\xi_x}^{GM}} $$

over the  range of integration for $k_x$. 


## Turbulent subrange

Instead, these observed turbulent subranges  could belong to a regime of layered anisotropic stratified turbulence \citep[LAST;][]{Riley2008, Falder2016}. Since there must be continuity between  LAST and  inertial-convective regimes, the dissipation rate of turbulent kinetic energy, $\epsilon$, can be estimated using the inertial-convective parametrisation for a passive scalar \citep{Sreenivasan1996}. If turbulent motions affect temperature in the same way as they affect density, this parametrization can be expressed as

$$ \left\langle \phi_{\xi_x} \right\rangle_p = \frac{4\pi\Gamma}{N^2} C_T \epsilon^{2/3} (2\pi k_x)^{+1/3}, $$

where $\Gamma$ is the turbulent flux coefficient and $C_T$ is the Obukhov-Corrsin constant \citep{Klymak2007turbulence, Falder2016}. Diapycnal diffusivity, $K$, is estimated from $\epsilon$ using the Osborn relationship

$$ K =  \frac{\frac{g}{\rho}\langle w'\rho' \rangle_p}{\frac{-g}{\rho}\frac{d\rho}{dz}} = \frac{B}{N^2} = \frac{\Gamma \epsilon}{N^2}, $$

 where $g$ is gravitational  acceleration, $\rho$ is  mean density, $\rho'$ is  density perturbation, $w'$ is vertical component of the velocity perturbation and $B$ is buoyancy flux \citep{Osborn1980}. Angular brackets denote averaging over time. This  equation  assumes a balance between turbulent production, dissipation, and work done against buoyancy for the case of a statistically steady flow. The Osborn relationship is also assumed in derivation of Equation (\ref{eq:tranche7and8_epsilon_from_turbulence}).

By combining  Equations (\ref{eq:tranche7and8_epsilon_from_turbulence}) and (\ref{eq:osborn_relation}), we obtain

$$ \log_{10}(K) = \frac{3}{2} \hat{c}  - \frac{1}{2} \log_{10}(\Gamma) + \log_{10}(N) - \frac{3}{2} \log_{10}(C_T) - \frac{7}{2}\log_{10}(2) - 2\log_{10}(\pi), $$

 where $\hat{c}$ is the intercept of a straight line of gradient $+1/3$ fitted in $\log_{10}(k_x)$&ndash;$\log_{10}\left\langle \phi_{\xi_x} \right\rangle_p$ space using least squares linear regression.


 Numerical uncertainties in estimated values of $K$ depend on uncertainties in $\hat{c}$, $\left\langle N \right\rangle_p$, $\Gamma$ and $C_T$. Error in the fitted intercept of a straight line of fixed gradient is given by
 
$$ \sigma_{\hat{c}} = \sqrt{\frac{1}{\sum_{i=1}^{D} \frac{1}{\sigma_i^2} }}, $$

where $\sigma_{\hat{c}}$ is the error on $\hat{c}$, $\sigma_i$ is the standard deviation on the $i$\textsuperscript{th} point of the mean spectrum $\left\langle \phi_{\xi_x} \right\rangle_p$, and $D$ is the total number of points in the mean spectrum. Estimated $\sigma_{\hat{c}}$ are in the range 0.017--0.1. 
Errors in $\hat{c}$ and $\left\langle N \right\rangle_p$ are combined in quadrature to yield

$$ (\sigma_K)^2 = \left(\frac{3}{2}\sigma_{\hat{c}}\right)^2 + \left(\log_{10}(e)\frac{\sigma \left\langle N \right\rangle_p }{\left\langle N \right\rangle_p}\right)^2, $$

where $\sigma_K$ and $\sigma \left\langle N \right\rangle_p$ are the uncertainties in $\log_{10}(K)$ and $\left\langle N \right\rangle_p$, respectively, and it has been assumed that $\sigma_{\hat{c}}$ and $\sigma \left\langle N \right\rangle_p $ are normally distributed. For $\sigma_{\hat{c}}$ in the range 0.017--0.1 and $\left\langle N \right\rangle_p = 2.0 \pm 0.3$~cph, errors on $\log_{10}(K)$ are in the range $\sim 0.07$--0.17 logarithmic units. If, instead of the value of $\left\langle N \right\rangle_p = 2.0 \pm 0.3$~cph estimated from CTD casts, a value of $\left\langle N \right\rangle_p = 2 \pm 1$~cph is used, uncertainties in $\log_{10}(K)$ are in the range $\sim 0.22$--0.26 logarithmic units.

In assessing the the effect of uncertainties in $\Gamma$ and $C_T$, we do not formally propagate errors because the underlying probability distributions are not accurately determined. Instead, the effects of choosing the lower and upper bounds given by \citet{Sreenivasan1996} and \citet{Mashayek2017} are explored.  Taking lower and upper bounds of $\Gamma = 0.1$ and $\Gamma = 0.4$, respectively, implies that assumption of a constant value of $\Gamma = 0.2$ introduces an uncertainty in $\log_{10}(K)$ of $\lesssim 0.15$ logarithmic units as shown by

$$ \log_{10}(K) = -\frac{1}{2}\log_{10}(\Gamma) +  A \\
 = \left( 0.5 + A \right) \ \mathrm{logarithmic \ units \ for} \ \Gamma = 0.1 \\
 \approx \left( 0.35 + A \right) \ \mathrm{logarithmic \ units \ for} \ \Gamma = 0.2 \\
 \approx \left( 0.2 + A \right) \ \mathrm{logarithmic \ units \ for} \ \Gamma = 0.4, $$
 
where $A$ includes all terms that are independent of $\Gamma$.
A similar analysis shows that, for $C_T$ in the range 0.3--0.5, assumption of constant $C_T = 0.4$ introduces uncertainties in $\log_{10}(K)$ of up to 0.1 logarithmic units.

## Notation

 Symbol | Description | Value | Unit | Dimension | Variable name
---|---|---|---|---|---
 $f$ | Coriolis parameter | | $\mathrm{s}^{-1}$ | $\mathrm{T}^{-1}$ | `f`
 $g$ | Gravitational acceleration | 9.81 | $\mathrm{m}$ $\mathrm{s}^{-2}$ | $\mathrm{L}$ $\mathrm{T}^{-1}$ | `g`
 $\hat{c}$ | Intercept fitted to identified turbulent spectral subranges |  |  | Dimensionless 
 $C_T$ | Obukhov-Corrsin constant | 0.4 |  | Dimensionless | `CT`
 $j_{res}$ | Parameter for multi-taper Fourier transforms |  |  | | `jres`
 $K$ | Diapycnal (eddy) diffusivity |  | $\mathrm{m}^2$ $\mathrm{s}^{-1}$ | $\mathrm{L}^2$ $\mathrm{T}^{-1}$ | `K`
  $K_{res}$ | Parameter for multi-taper Fourier transforms |  |  | | `Kres`
 $k_x$ | Component of horizontal wavenumber in plane of seismic image |  | cpm |$\mathrm{L}^{-1}$ | `kx`
 $\log_{10}(k_x)$ |  |  |  |  | `log_kx`
 $\log_{10}(\phi_{\xi_x})$ |  |  |  |  | `log_phix`
 $\left\langle \phi_{\xi_x} \right\rangle$ |  |  |  |  | `mean_phix`
 $\sigma\left\langle \phi_{\xi_x} \right\rangle$ |  |  |  |  | `std_phix`
 $\log_{10}\left\langle\phi_{\xi_x}\right\rangle$ |  |  |  |  | `log_mean_phix`
 $\log_{10}(e) \frac{\sigma\left\langle\phi_{\xi_x}\right\rangle}{\left\langle\phi_{\xi_x}\right\rangle}$ |  |  |  |  | `std_log_mean_phix`
 $\left\langle\log_{10}(\phi_{\xi_x})\right\rangle$ |  |  |  |  | `mean_log_phix`
 $\sigma\left\langle\log_{10}(\phi_{\xi_x})\right\rangle$ |  |  |  |  | `std_log_phix`
 $\Gamma$ | Turbulent flux coefficient | 0.2 |  | Dimensionless | `gamma`
 $\epsilon$ | Dissipation rate of turbulent kinetic energy |  | $\mathrm{m}^2$ $\mathrm{s}^{-3}$ | $\mathrm{L}^2$ $\mathrm{T}^{-3}$ | `epsilon`
 $\phi_\xi(k_x)$ | Power spectrum of vertical displacement |  | $\mathrm{m}^2$ $[\mathrm{cpm}]^{-1}$ | $\mathrm{L}^{3}$ | `phi`
 $\phi_{\xi_x}(k_x)$ | Power spectrum of horizontal gradient of vertical displacement |  | $[\mathrm{cpm}]^{-1}$ | $\mathrm{L}$ | `phix`
 $\left\langle \phi_{\xi_x} \right\rangle_p$ | Mean $\phi_{\xi_x}$ within seismically imaged pycnocline |  | $[\mathrm{cpm}]^{-1}$ | $\mathrm{L}$ 
 $\phi_{\xi_x}^{\mathrm{DDT}}(k_x)$ | Direct data transform | | Arbitrary units 
 $\phi_{\xi_x}^{GM}(k_x)$ | Garrett-Munk spectrum of horizontal gradient of vertical displacement | | [cpm]$^{-1}$ | $\mathrm{L}$ 
