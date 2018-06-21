from setuptools import setup

setup(
    name='g4netx',
    version='0.1',
    description=(
        'Find all overlapping G4 patterns using network analysis'
    ),
    author='Matthew Parker',
    packages=[
        'g4netx',
    ],
    install_requires=[
        'networkx',
    ],
)