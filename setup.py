from setuptools import find_packages, setup
from version import version as this_version


try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

import os
import sys

# Utility function to read the README file.
# Used for the long_description.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# parse_requirements() returns generator of pip.req.InstallRequirement objects
install_reqs = parse_requirements('requirements.txt', session='hack')

reqs = [str(ir.req) for ir in install_reqs]

long_description = "!!!!! missing pandoc do not upload to PyPI !!!!"
try:
    import pypandoc
    # TODO needs pandoc installed
    long_description = pypandoc.convert('README.md', 'rst')
except ImportError:
    print("Could not import pypandoc - required to package epitome", file=sys.stderr)
except OSError:
    print("Could not convert - pandoc is not installed", file=sys.stderr)

setup(
    name='epitome',
    version=this_version,
    description='epigenetic learning',
    author='Alyssa Morrow',
    author_email='akmorrow@berkeley.edu',
    url="https://github.com/akmorrow13/epitome",
    install_requires=reqs,
    long_description=long_description,
    packages=find_packages(exclude=['*.test.*']))