from setuptools import find_packages, setup

with open("README.md", "r") as readme:
    long_description = readme.read()

setup(
    name = 'uwTOPSIS',
    packages = find_packages(include=['uwTOPSIS']),
    version = '0.1.0',
    author = 'Aaron Lopez-Garcia',
    author_email='aaron.lopez@uv.es',
    description = 'Unweighted TOPSIS method',
    long_description=long_description,
    license = 'MIT',
    install_requires=['pandas >= 1.2.4',
                      'numpy >= 1.19',
                      'scipy >= 1.6.3'],
    classifiers=["Programming Language :: Python :: 3.8.5"],
)
