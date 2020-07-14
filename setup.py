import subprocess
import sys
from distutils.errors import DistutilsOptionError
from distutils.spawn import find_executable
from os import path
from pathlib import Path

from setuptools import find_packages, setup
from setuptools.command.develop import develop
from setuptools.command.install import install

# Speed up entrypoints
# see https://github.com/ninjaaron/fast-entry_points
# import fastentrypoints

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
    "xarray>=0.16.0",
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
    "napari==0.3.1",
    "tqdm>=4.46",
    "click>=7.1.2,<8",
    "PyYAML>=5.3.1",
    "sphinx>=3.0.3",
    "sphinx_rtd_theme>=0.4.3",
    "sphinx-autodoc-typehints>=1.10.3",
    "numpydoc>=0.9.2",
    "requests",
    "jupyter",
    "pylint",
    "isort",
    "black",
    "strictyaml"
]

test_requirements = ["pytest", "pytest-xdist", "pytest-cov", "pytest-datadir"]


class InstallMatlabEngineMixin:
    cmd_options = [("matlab-root=", None, "MATLAB installation directory (MATLABROOT)")]

    def _initialize(self):
        self.matlab_root = None

    def _finalize(self):
        if self.matlab_root is not None:
            self.matlab_root = Path(self.matlab_root)
            if not self.matlab_root.is_dir():
                raise DistutilsOptionError(
                    'MATLAB installation directory "{}" does not exist'.format(
                        self.matlab_root
                    )
                )

    def _install_matlab_engine(self):
        if self.matlab_root:
            matlab_root = self.matlab_root
            install_failed_error = True
        else:
            matlab_exe = find_executable("matlab")
            if not matlab_exe:
                return

            matlab_exe_parent = Path(matlab_exe).parent
            install_failed_error = False

            if matlab_exe_parent.name == "bin":
                matlab_root = matlab_exe_parent.parent
            else:
                # /matlabroot/bin/win64
                matlab_root = matlab_exe_parent.parent.parent

            print('MATLAB was found. MATLABROOT: "{}"'.format(matlab_root))

        engine_dir = matlab_root / "extern" / "engines" / "python"

        if not engine_dir.is_dir():
            raise EnvironmentError(
                'Cannot find "MATLAB engine for Python" in MATLAB root "{}"'.format(
                    matlab_root
                )
            )

        print('Installing "MATLAB engine for Python" from "{}"...'.format(engine_dir))

        engine_setup_path = str(engine_dir / "setup.py")
        install_command = [sys.executable, engine_setup_path, "install"]

        p = subprocess.run(
            install_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        if p.returncode != 0:
            try:
                raise EnvironmentError(
                    'An error occurred while installing "MATLAB engine for Python"\n'
                    + "{}\n".format(p.stderr.decode("utf-8"))
                )
            except EnvironmentError as err:
                if install_failed_error:
                    raise
                else:
                    print(err, file=sys.stderr)
        else:
            print(
                '{}\n\n"Matlab Engine for Python" is successfully installed'.format(
                    p.stdout.decode("utf-8")
                )
            )


class InstallCommand(install, InstallMatlabEngineMixin):
    user_options = install.user_options + InstallMatlabEngineMixin.cmd_options

    def initialize_options(self):
        install.initialize_options(self)
        self._initialize()

    def finalize_options(self):
        install.finalize_options(self)
        self._finalize()

    def run(self):
        self._install_matlab_engine()
        install.run(self)


class DevelopCommand(develop, InstallMatlabEngineMixin):
    user_options = develop.user_options + InstallMatlabEngineMixin.cmd_options

    def initialize_options(self):
        develop.initialize_options(self)
        self._initialize()

    def finalize_options(self):
        develop.finalize_options(self)
        self._finalize()

    def run(self):
        self._install_matlab_engine()
        develop.run(self)


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
    packages=find_packages(include=["pharedox", "pharedox.*"]),
    python_requires=f">={MIN_PY_VER}",
    install_requires=requirements,
    project_urls={"Bug Reports": "https://github.com/ApfeldLab/pharynx_redox/issues"},
    entry_points="""
        [console_scripts]
        pharedox=pharedox.scripts.__main__:cli
    """,
    cmdclass={"install": InstallCommand, "develop": DevelopCommand},
)
