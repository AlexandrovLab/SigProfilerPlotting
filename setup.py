from setuptools import setup

def readme():
	with open('README.rst') as f:
		return(f.read())


setup(name='sigProfilerPlotting',
		version='0.1.5',
		description='SigProfiler plotting tool',
		url='',
		author='Erik Bergstrom',
		author_email='ebergstr@eng.ucsd.edu',
		license='UCSD',
		packages=['sigProfilerPlotting'],
		install_requires =[
			"matplotlib"],
		include_package_data=True,
		zip_safe=False)