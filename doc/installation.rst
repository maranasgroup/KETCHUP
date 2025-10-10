Command line installation instructions for KETCHUP
==================================================

Command prompt lines are set within <>. File names are set within {}. Directories are set within “”. When sending command, do not type in <, >, {, or }.

Requirements for KETCHUP:

* Python 3.12.9
* Pyomo 6.9.4
* Pandas 2.3.3
* Openpyxl 3.1.5
* Scipy 1.16.2
* Cobra 0.29.1
* Pyyaml 6.0.3

Requirements for manual compilation of IPOPT with HSL solver
============================================================

* LAPACK 3.11.0
* HSL -Galahad package 5.10000.0
* Ipopt 3.14.19
* ThirdParty-ASL-stable-2.0

Installation of KETCHUP Environment for Anaconda
================================================

For convenience an `environment YAML file <../installation/ketchup_environment.yml>`__  is provided in the installation directory.
To use it, download the file and then run::


   conda env create -f ketchup_environment.yml


which will create the pyomo_ketchup anaconda environment. You might need to add in the path to where you downloaded the environment YAML file.

Command line Pyomo with IPOPT from scratch.
===========================================


#. Create conda environment

   ``<conda create –-name pyomo_ketchup python=3.12.9>``

#. Activate environment

   ``<conda activate pyomo_ketchup>``

#. Install pandas 2.3.3

   * ``<conda install pandas=2.3.3>``

   or

   * ``<pip install pandas==2.3.3>``

#. Install openpyxl 3.1.5

   * ``<conda install openpyxl=3.1.5>``

   or

   * ``<pip install openpyxl==3.1.5>``

#. Install IPOPT 3.14.19

   * ``<conda install -c conda-forge ipopt=3.14.19>``

   or

   * ``<pip install ipopt==3.14.19>``

#. Install Pyomo 6.9.4

   * ``<conda install -c conda-forge pyomo=6.9.4>``

   or

   * ``<pip install pyomo==6.9.4>``

#. Install scipy 1.16.2

   * ``<conda install scipy=1.16.2>``

   or

   * ``<pip install scipy==1.16.2>``

#. Install pyyaml 6.0.3

   * ``<conda install -c conda-forge pyyaml=6.0.3>``

   or

   * ``<pip install pyyaml==6.0.3>``

#. Install cobra 0.29.1

   ``<pip install cobra==0.29.1>``

#. Optional: Install streamlit 1.50.0 (for GUI)

   ``<pip install streamlit==1.50.0>``

\*OPTIONAL\* installation of AMPL
=================================

AMPL module is to allow for pyomo to interface with alternative solvers such as BONMIN (MINLP) or CONOPT(NLP, license required); KETCHUP will still run without this module installed. We have tested this and confirmed that AMPL currently has lower performance with Ipopt HSL solvers than manually compiling Ipopt with HSL libraries.

#. Install amplpy 0.14.0

   * ``<conda install amplpy=0.14.0>``

   or

   * ``<pip install amplpy==0.14.0>``

#. Install IPOPT module with amply (license required, academic license provides free solvers, see screenshots below for ampl)

   #. Install amplpy modules

      * ``<python -m amplpy.modules install coin>``

   #. Activate amplpy license

      * ``<python -m amplpy.modules activate LICENSE>``

Linking HSL Libraries with IPOPT
================================

Choose a directory and set it to $ENV, this will be the directory used for compiling HSL solver files. This installation assumes that you are already working in the anaconda environment you created above.

#. Apply for an academic license for HSL package at <https://www.hsl.rl.ac.uk/>

   * Download file {hsl-galahadd-5.10000.0.tar.gz}

#. Install meson 1.7.2 for HSL

   * ``<pip install meson==1.7.2>``

#. Change to $ENV directory

   * ``<cd $ENV>``

#. Extract HSL files (there will be a resulting folder called “hsl\_subset” upon extraction)

   * ``<tar xf hsl-galahadd-5.10000.0.tar.gz>``

#. Change directory into the “hsl\_subset” directory

   * ``<cd hsl\_subset>``

#. Set environment variable HSLSUBSET to the current directory (find your absolute directory path with <pwd> command)

   * ``<export HSLSUBSET=absolute-path-to-current-directory>``

#. In the “makefiles” directory, modify script “pc64.lnx.gfo” to select your proper gfortran and c compiler. See screenshots below.

   * ![](data:image/png;base64...)

#. In current directory, modify meson\_options.txt file to turn off modules option

   * ![](data:image/png;base64...)

#. Move to the “src” directory and build the fortran libraries for both single and double precisions

   #. ``<cd src>``
   #. ``<make -s -f $HSLSUBSET/makefiles/pc64.lnx.gfo hsl PRECIS=single>``
   #. ``<make -s -f $HSLSUBSET/makefiles/pc64.lnx.gfo hsl PRECIS=double>``

#. Setup meson build in main directory

   #. ``<cd $HSLSUBSET>``
   #. ``<meson setup builddir>``

#. Move all mod files from “modules” directory into “builddir”

   #. ``<cp $HSLSUBSET/modules/pc64.lnx.gfo/single/hsl\_\* $HSLSUBSET/builddir>``
   #. ``<cp $HSLSUBSET/modules/pc64.lnx.gfo/double/hsl\_\* HSLSUBSET/builddir>``

#. Compile files

   * ``<meson compile -C builddir>``

#. Rename libhsl\_subset.so found in “builddir” to libhsl.so, and append the “builddir” directory to LD\_LIBRARY\_PATH

   #. ``<mv ./builddir/libhsl\_subset.so ./builddir/libhsl.so>``
   #. ``<export LD\_LIBRARY\_PATH=$LD\_LIBRARY\_PATH:/$HSLSUBSET/builddir>``

Manual compilation of Ipopt with HSL.

#. Go to working environment $ENV set when compiling HSL

   * ``<cd $ENV>``

#. Clone Ipopt from github

   * ``<git clone <https://github.com/coin-or/Ipopt.git>>``

#. Clone ThirdParty-ASL from github (this package is for interfacing IPOPT with Pyomo)

   * ``<git clone <https://github.com/coin-or-tools/ThirdParty-ASL> -b stable/2.0>``

#. Move to “ThirdParty-ASL” directory, configure, and compile

   #. ``<cd ThirdParty-ASL>``
   #. ``<./get.ASL>``
   #. ``<./configure -prefix=$ENV>``
   #. ``<make>``
   #. ``<make install>``

#. Move to Ipopt directory

   * ``<cd $ENV/Ipopt>``

#. Create a build directory to for compiling and creating of ipopt executable

   #. ``<mkdir build>``
   #. ``<cd build>``

#. Configure build with openmp flags. Openmp allows for parallel computing that speeds-up solve times for ipopt

   #. (if compiler is mkl may need to module load mkl first)

      * ``<../configure --prefix=$ENV ADD\_CFLAGS=-fopenmp ADD\_FFLAGS=-fopenmp ADD\_CXXFLAGS=-fopenmp>``

   #. Make,test, and build

      #. ``<make>``
      #. ``<make test>``
         *. With HSL solver in place, all test should pass
      #. ``<make install>``

#. Add PATH to the Ipopt executable

   #. ``<export PATH=$PATH:$ENV/Ipopt/build/src/Apps/AmplSolver/>``

Accessing AMPL python installation

#. Sign up for an account at https://portal.ampl.com/user/ampl/home
#. Login to AMPL, and select “Download or Modify your AMPL CE license”
#. Select desired solvers and select Python for installation information (this will include license information)
