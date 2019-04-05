from setuptools import setup
from setuptools.command.install import install
import warnings


def readme():
	with open('README.rst') as f:
		return(f.read())


#Set up the machinery to install custom fonts.  Subclass the setup tools install 
#class in order to run custom commands during installation.  
class move_ttf(install):
    def run(self):
        """
        Performs the usual install process and then copies the True Type fonts 
        that come with clearplot into matplotlib's True Type font directory, 
        and deletes the matplotlib fontList.cache 
        """
        #Perform the usual install process
        install.run(self)
        #Try to install custom fonts
        try:
            import os, shutil
            import matplotlib
            import matplotlib.font_manager

            #Find where matplotlib stores its True Type fonts
            mpl_data_dir = os.path.dirname(matplotlib.matplotlib_fname())
            mpl_ttf_dir = os.path.join(mpl_data_dir, 'fonts', 'ttf')
            cp_ttf_dir = os.path.dirname(os.path.realpath(__file__))

            file_names = ["Times New Roman.ttf", "Arial.ttf", "Courier New.ttf", 
            			  "Courier New Bold.ttf", "Arial Bold.ttf", "Times New Roman Bold.ttf"]
            for file in file_names:
            	old_path = os.path.join(cp_ttf_dir, "fonts/" + file)
            	new_path = os.path.join(mpl_ttf_dir, file)
            	shutil.copyfile(old_path, new_path)
            matplotlib.font_manager._rebuild()

        except:
            warnings.warn("WARNING: An issue occured while installing the fonts.")

        




setup(name='sigProfilerPlotting',
		version='0.1.26',
		description='SigProfiler plotting tool',
		url='',
		author='Erik Bergstrom',
		author_email='ebergstr@eng.ucsd.edu',
		license='UCSD',
		packages=['sigProfilerPlotting'],
		install_requires =[
			"matplotlib"],
		include_package_data=True,
	    #Specify the custom install class
	    cmdclass={'install' : move_ttf},
		zip_safe=False)

