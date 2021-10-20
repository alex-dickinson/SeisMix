<!---[![Chat on Gitter](https://img.shields.io/gitter/room/calliope-project/calliope.svg?style=flat-square)](https://gitter.im/calliope-project/calliope)--->

[![Master branch build status](https://img.shields.io/azure-devops/build/calliope-project/371cbbaa-fa6b-4efb-9b23-c4283a8e33eb/1?style=flat-square)](https://dev.azure.com/calliope-project/calliope/_build?definitionId=1)
[![Documentation build status](https://img.shields.io/readthedocs/calliope.svg?style=flat-square)](https://readthedocs.org/projects/calliope/builds/)
[![Test coverage](https://img.shields.io/codecov/c/github/calliope-project/calliope?style=flat-square&token=b4fd170f0e7b43679a8bf649719e1cea)](https://codecov.io/gh/calliope-project/calliope)
[![PyPI version](https://img.shields.io/pypi/v/calliope.svg?style=flat-square)](https://pypi.python.org/pypi/calliope)
[![Anaconda.org/conda-forge version](https://img.shields.io/conda/vn/conda-forge/calliope.svg?style=flat-square&label=conda)](https://anaconda.org/conda-forge/calliope)
[![JOSS DOI](https://img.shields.io/badge/JOSS-10.21105/joss.00825-green.svg?style=flat-square)](https://doi.org/10.21105/joss.00825)


# Mixing codes for Kathy

Make sure code works with example data from FSC

Unused modules?:
contours_to_midpoints_fs3d.py
contours_to_midpoints.py


install with pip:
Make venv
pip install .

To install optional packages for running example jupyter notebooks:
```bash 
pip install .[examples]
```

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

### Using pip

Set up virtual pip environment.

```bash
pip install SeisMix -r requirements.txt
```



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




