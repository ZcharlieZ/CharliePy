import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "CharliePy",
    version = "1.0.0",
    author = "Carlo Cannarozzo",
    author_email = "carlo.cannarozzo3@unibo.it",
    description = ("CharliePy modules"),
    #url = "http://packages.python.org/an_example_pypi_project",
    packages=["Analysis_Catalogue"],
    long_description=read('README'),
)