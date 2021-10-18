import setuptools
import versioneer

with open("README.md", "r") as fh:
    long_description = fh.read()
with open("requirements.txt", "r") as fh:
    requirements = [line.strip() for line in fh]

setuptools.setup(
    name="SeisMix",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Alex Dickinson",
    author_email="nad38@cantab.ac.uk",
    description="WHAT IS THE REASON YOU BUILD THIS PROJECT AND WHAT IT DOES",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url = "https://github.com/alex-dickinson/SeisMix",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=requirements,
)