from setuptools import setup, find_packages
from pathlib import Path

here = Path(__file__).resolve().parent
README = (here / "README.md").read_text(encoding="utf-8")
VERSION = (here / 'hiergo' / "VERSION").read_text(encoding="utf-8").strip()

setup(
    name='hiergo',
    packages=['hiergo',
              ] + find_packages(exclude=['tests', 'tests.*']),
    include_package_data=True,
    entry_points={
        "console_scripts": ["hiergo=hiergo.cli:execute_cli"],
    },
    version=VERSION,
    license='gpl-3',
    description='Structural atomistic models of graphene oxide via a modular tiling approach',
    long_description=README,
    long_description_content_type='text/markdown',
    author='Natalya A Garcia, Joel B Awuah, Chaoyue Zhao, Filip Vuković, Tiffany R Walsh',
    author_email='sherif.tawfic@gmail.com',
    url='https://github.com/sheriftawfikabbas/joltml',
    keywords=['ai', 'machine learning', 'machine learning workflow'
              'model tracking'],
    install_requires=['xgboost',
                      'pandas',
                      'numpy',
                      'torch',
                      'tensorflow',
                      'sklearn'],

)