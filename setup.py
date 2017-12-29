import os
from setuptools import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()


setup(name='mixcurve',
      version='0.1',
      description='Package for interpreting 1:1 Mixing Studies',
      url='',
      author='Nik Krumm',
      author_email='nkrumm@gmail.com',
      license='MIT',
      packages=['mixcurve'],
      include_package_data=True,
      zip_safe=False,
      install_requires=required)