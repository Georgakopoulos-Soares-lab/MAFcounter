from setuptools import setup, find_packages

setup(
    name="maf_counter",
    version='0.2',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[],
    entry_points={
        "console_scripts": [
            "maf_counter=maf_counter.cli:main",
        ],
    },
    package_data={
        "maf_counter": ["bin/maf_counter"],  # Include the binary in the package
    },
     author='Patsakis Michail , Provatas Kimon, Ilias Georgakopoulos Soares, Ioannis Mouratidis',  # Your name
    author_email='kap6605@psu.edu , mpp5977@psu.edu',  # Your email
    description="A high-performance k-mer counting tool for MAF alignments.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Georgakopoulos-Soares-lab/MAFcounter",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
)
