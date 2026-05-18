from pathlib import Path

from setuptools import find_packages, setup

README = Path("README.md")

setup(
    name="lipidbuilder",
    version="0.1.0",
    author="Julianne Hoeflich",
    description="Build and parameterize lipid bilayer systems with Packmol and OpenFF.",
    long_description=README.read_text() if README.exists() else "",
    long_description_content_type="text/markdown",
    packages=find_packages(exclude=["lipidbuilder.results*", "lipidbuilder.tests*"]),
    include_package_data=True,
    package_data={
        "lipidbuilder": [
            "data/available-molecules/PulledLipid.csv",
            "data/available-molecules/lipid_pdbs/*/*.pdb",
            "data/available-molecules/solvent_pdbs/*.pdb",
            "data/forcefields/*.offxml",
            "data/gromacs/*.mdp",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10,<3.12",
)
