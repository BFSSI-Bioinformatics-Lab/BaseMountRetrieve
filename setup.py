import setuptools
from BaseMountRetrieve.basemountretrieve import __version__, __author__, __email__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="BaseMountRetrieve",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="Internal BFSSI package for retrieving files from BaseMount",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bfssi-forest-dussault/BaseMountRetrieve",
    packages=setuptools.find_packages(),
    scripts=['BaseMountRetrieve/basemountretrieve.py'],
    install_requires=['click', 'pandas', 'numpy']
)
