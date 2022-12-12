from setuptools import setup, find_packages

VERSION = '1.3'
DESCRIPTION = 'Redshift a radio source to a higher redshift'

with open('README.md') as f:
   LONG_DESCRIPTION = f.read()

# Setting up
setup(
    # the name must match the folder name 'verysimplemodule'
    name="redshifting",
    version=VERSION,
    author="Jurjen de Jong",
    author_email="<jurjong@proton.me>",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=['astropy', 'reproject', 'matplotlib', 'numpy'],
    keywords=['python', 'redshifting', 'radio astronomy'],
    classifiers=[]
)