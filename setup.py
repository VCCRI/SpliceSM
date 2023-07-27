from setuptools import setup, find_packages
import io

setup(
    name='Sgen',
    description='Predicting splice-altering variants',
    version='0.1.0',
    author='Steve Monger',
    author_email='s.monger@victorchang.edu.au',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    #packages=['sgen'],
    scripts=['bin/run.py'],
    #url='http://pypi.python.org/pypi/temp',
    license='MIT',
    install_requires=[
        "BioPython",
        "Cython",
    ],
)
