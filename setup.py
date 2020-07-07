from setuptools import setup, find_packages
from os import path
import sys

# Speed up entrypoints
# see https://github.com/ninjaaron/fast-entry_points
try:
    import fastentrypoints
except (ImportError, ModuleNotFoundError):
    from setuptools.command import easy_install
    import pkg_resources
    easy_install.main(['fastentrypoints'])
    pkg_resources.require('fastentrypoints')
    import fastentrypoints

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

MIN_PY_MAJOR_VER = 3
MIN_PY_MINOR_VER = 7
MIN_PY_VER = f"{MIN_PY_MAJOR_VER}.{MIN_PY_MINOR_VER}"
DISTNAME = "pharedox"
DESCRIPTION = "Software tools for measuring cytosolic redox state in the pharynx of the nematode C. elegans"
LONG_DESCRIPTION = long_description
LICENSE = "BSD 3-Clause"
DOWNLOAD_URL = "https://github.com/ApfeldLab/pharynx_redox/"

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Environment :: X11 Applications :: Qt",
    "Intended Audience :: Education",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Visualization",
    "Topic :: Scientific/Engineering :: Information Analysis",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Utilities",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: MacOS",
]

if sys.version_info < (MIN_PY_MAJOR_VER, MIN_PY_MINOR_VER):
    sys.stderr.write(
        f"You are using Python "
        f"{'.'.join(str(v) for v in sys.version_info[:3])}.\n\n"
        f"napari only supports Python {MIN_PY_VER} and above.\n\n"
        f"Please install Python {MIN_PY_VER} or later.\n"
    )
    sys.exit(1)

requirements = [
    "numpy>=1.18",
    "xarray==0.15.1",
    "netCDF4==1.5.3",
    "scikit-image==0.17.2",
    "tifffile==2020.5.11",
    "scipy>=1.4.1",
    "numba>=0.49.1",
    "scikit-learn>=0.23.0",
    "pandas>=1.0.3,<2",
    "statsmodels==0.11.1",
    "simpleITK==1.2.4",
    "matplotlib>=3.2.1",
    "seaborn>=0.10.1",
    "PyQt5>=5.14.2",
    "napari>=0.3.1",
    "tqdm>=4.46",
    "click>=7.1.2,<8",
    "strictyaml>=1.0.6",
    "sphinx>=3.0.3",
    "sphinx_rtd_theme>=0.4.3",
    "sphinx-autodoc-typehints>=1.10.3",
    "numpydoc>=0.9.2",
    "requests>2",
    "isort<5",
    "jupyter",
    "pylint",
    "black",
]

test_requirements = ["pytest", "pytest-xdist", "pytest-cov", "pytest-datadir"]

setup(
    name="pharedox",
    version="0.0.1",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=DOWNLOAD_URL,
    author="Sean Johnsen",
    author_email="sean.b.johnsen@gmail.com",
    maintainer="Sean Johnsen",
    maintainer_email="sean.b.johnsen@gmail.com",
    classifiers=CLASSIFIERS,
    keywords="image-analysis biology microscopy",
    packages=find_packages(),
    python_requires=f">={MIN_PY_VER}",
    install_requires=requirements,
    setup_requires=["pytest-runner",],
    tests_require=test_requirements,
    test_suite="tests",
    project_urls={"Bug Reports": "https://github.com/ApfeldLab/pharynx_redox/issues"},
    entry_points="""
        [console_scripts]
        pharedox=pharedox.scripts.__main__:cli
    """,
)
