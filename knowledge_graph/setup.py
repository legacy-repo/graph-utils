import os
from setuptools import setup

version = "1.0.0"

with open("README.md", "r") as fh:
    long_description = fh.read()


with open("requirements.txt", "r") as f:
    requirements = f.readlines()


def get_packages(package):
    """Return root package and all sub-packages."""
    return [dirpath
            for dirpath, dirnames, filenames in os.walk(package)
            if os.path.exists(os.path.join(dirpath, "__init__.py"))]


setup(
    name="knowledge-graph",
    version=version,
    author="Jingcheng Yang",
    author_email="yjcyxky@163.com",
    description="A Python project that allows you to analyse proteomics and clinical data, and integrate and mine knowledge from multiple biomedical databases widely used nowadays.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/biomedgps/biomedgps.git",
    packages=get_packages("builder") + get_packages("importer"),
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'knowledge-graph=builder:knowledge_graph',
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='==3.7.9',
    include_package_data=True,
)
