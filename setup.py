from setuptools import setup

def readme():
	with open('README.rst') as f:
		return(f.read())


setup(name='sigProfilerPlotting',
		version='0.1',
		description='SigProfiler plotting tool',
		url='',
		author='Erik Bergstrom/Maria',
		author_email='ebergstr@eng.ucsd.edu',
		license='UCSD',
		packages=['sigProfilerPlotting'],
		zip_safe=False)