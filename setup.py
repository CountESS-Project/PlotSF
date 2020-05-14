import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = ["matplotlib >= 3.0.0", "pandas >= 1.0.0", "numpy"]

setuptools.setup(
    name="PlotSF",
    version="0.0.1",
    author="Alan F Rubin",
    author_email="alan.rubin@wehi.edu.au",
    description="A plotting library for multiplexed assay of variant effect data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CountESS-Project/PlotSF",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=requirements,
)
