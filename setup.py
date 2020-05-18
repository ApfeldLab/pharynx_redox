from setuptools import setup, find_packages
from os import path
import sys

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
DOWNLOAD_URL = "https://github.com/ApfeldLab/pharedox/"

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Environment :: X11 Applications :: Qt",
    "Intended Audience :: Education",
    "Intended Audience :: Science/Research",
    # 'License :: OSI Approved :: BSD License',
    # 'Programming Language :: C',
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
    "numpy",
    "xarray",
    "netcdf4",
    "scikit-image",
    "scipy",
    "numba",
    "scikit-learn",
    "pandas",
    "statsmodels",
    "simpleITK",
    "matplotlib",
    "seaborn",
    "pyqt5",
    "pyqtgraph",
    "napari",
    "tqdm",
    "pyyaml",
    "sphinx",
    "sphinx_rtd_theme",
    "sphinx-autodoc-typehints",
    "numpydoc",
    "cached-property",
    "jupyter",
    "pylint",
    "black",
]

test_requirements = [
    "pytest",
    "pytest-xdist",
    "pytest-cov",
]

setup(
    name="pharedox",
    version="0.0.1",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    download_url=DOWNLOAD_URL,
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
    project_urls={"Bug Reports": "https://github.com/ApfeldLab/pharedox/issues"},
)
