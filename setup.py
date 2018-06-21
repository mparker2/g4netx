from setuptools import setup


setup(
    name='g4netx',
    version='0.1',
    description=(
        'Find all overlapping G4 patterns using network analysis'
    ),
    author='Matthew Parker',
    entry_points={
        'console_scripts': [
            'g4netx = g4netx.g4netx:g4netx_cli'
        ]
    },
    packages=[
        'g4netx',
    ],
    install_requires=[
        'networkx',
    ],
)