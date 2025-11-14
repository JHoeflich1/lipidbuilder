from setuptools import setup, find_packages
import os
setup(
    name="lipidbuilder",
    version="0.1.0",
    author="Julianne Hoeflich",
    description="A toolkit for building and parameterizing lipid bilayers with openFF.",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    url="https://github.com/JHoeflich1/lipidbuilder", 
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)