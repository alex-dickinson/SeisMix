# Mixing codes for Kathy

Make sure code works with example data from FSC

Unused modules?:
contours_to_midpoints_fs3d.py
contours_to_midpoints.py


install with pip:
Make venv
pip install -r requirements.txt

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

gsw?


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
- [`jupyter`](https://jupyter.org/)
- [`matplotlib`](https://matplotlib.org/)



