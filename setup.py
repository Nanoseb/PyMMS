import setuptools
import PyMMS

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyMMS",
    version=PyMMS.__version__,
    author="Sebastien Lemaire",
    author_email="sebastien.lemaire@soton.ac.uk",
    description="Generation of Unsteady RANS manufactured solution for CFD",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nanoseb/PyMMS",
    packages=setuptools.find_packages(),
    requires=['Sympy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires='>=3.2',
)
