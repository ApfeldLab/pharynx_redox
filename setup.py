from setuptools import setup
from os import path, system

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(
    name="pharynx-redox",
    version="0.1.0",
    description="Software tools for measuring cytosolic redox state in the pharynx of the nematode C. elegans",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/half-adder/pharynx_redox",
    author="Sean Johnsen",
    author_email="sean.b.johnsen@gmail.com",
    maintainer="Sean Johnsen",
    maintainer_email="sean.b.johnsen@gmail.com",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.7",
    ],
    keywords="celegans nematode redox image-analysis biology microscopy",
    packages=["pharynx_redox"],
    python_requires="~=3.7",
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "jupyter",
        "xarray",
        "tqdm",
        "statsmodels",
        "matplotlib",
        "xlrd",
        "elasticdeform",
        "attrs",
        "scikit-image",
        "scikit-learn",
        "seaborn",
        "scikit-fda @ git+https://github.com/half-adder/scikit-fda.git@develop#egg=scikit-fda",
        "tabulate",
        "cached-property",
    ],
    test_suite="tests",
    extras_require={
        "test": ["pytest", "pytest-cov"],
        "docs": ["sphinx", "sphinx-rtd-theme", "sphinx-autodoc-typehints"],
    },
    project_urls={"Bug Reports": "https://github.com/half-adder/pharynx_redox/issues"},
)
