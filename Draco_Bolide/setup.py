import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="DRACO-Bolide-Mucephie",
    version="0.0.1",
    author="June Parsons",
    author_email="mucephie@my.yorku.ca",
    description="Data Reduction at the Allan Carswell Observatory. This package contains tools that allow for the the automatic reduction and assisted analysis of CCD data ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Mucephie/DRACO",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)