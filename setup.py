"""SeisMix

"""

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()
with open("requirements.txt", "r") as fh:
    requirements = [line.strip() for line in fh]
with open("optional_requirements.txt", "r") as fh:
    optional_requirements = [line.strip() for line in fh]

setup(
    name="SeisMix",
    version="0.0.1.dev",
    author="Alex Dickinson",
    author_email="nad38@cantab.ac.uk",
    description="WHAT IS THE REASON YOU BUILD THIS PROJECT AND WHAT IT DOES",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/alex-dickinson/SeisMix",
    packages=find_packages(),
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific"
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
	extra_requirements=optional_requirements,
)
