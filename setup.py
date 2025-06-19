from setuptools import setup, find_packages
import os

setup(
    name='qchem',
    version=os.environ.get("PACKAGE_VERSION", "0.0.0"),
    packages=find_packages(),
    description='A Quantum Chemistry Package for Python. Used and Created by the UW Formula Nano Team to design a Nano Car for the 2026/7 Nano Car Race',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='UWFormulaNano',
    author_email='uwformulanano@gmail.com',
    url='https://github.com/UWFormulaN/qchem',
    install_requires=[
        'numpy>=1.26.4',  # Add your dependencies here
        'pandas>=2.2.2',

    ],
)