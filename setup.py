import os
import setuptools

name = 'stretchdyn'
version = '0.1.0'

with open('README.md', 'r') as fh:
    long_description = fh.read()

PACKAGE_PATH = os.path.abspath(os.path.join(__file__, os.pardir))

setuptools.setup(
    url='https://github.com/akey7/StretchDyn',
    name=name,
    version=version,
    author='Alicia Key',
    author_email='alicia@akey7.com',
    description='StretchDyn',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(PACKAGE_PATH, "test"),
    install_requires=[
        'pytest',
        'black',
        'mypy',
        'numpy',
    ]
)
