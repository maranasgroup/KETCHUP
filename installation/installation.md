# Command line for pyomo installation with IPOPT from scratch

Command prompt lines are set within &lt;&gt;. File names are set within {}. When sending command, do not type in &lt;, &gt;, {, or }.

# Requirements for simple installation
<ul>
<li>Python 3.8.5</li>
<li>Pyomo 6.4.0</li>
<li>Pandas 1.2.5</li>
<li>Openpyxl 3.0.7</li>
<li>Cobra 0.26.3</li>
<li>Scipy 1.8.0</li>
<li>Ipopt 3.14.0</li>
<li>Numpy 1.23.1</li>
</ul>

# Requirements for manual compilation of Ipopt with HSL solver
<ul>
<li>Coinhsl-2021.05.05.tar.gz</li>
<li>Ipopt-releases-3.14.0.tar.gz</li>
<li>ThirdParty-HSL-stable-2.1</li>
<li>ThirdParty-ASL-stable-2.0</li>
<li>HSL package 2021.05.05</li>
<li>LAPACK 3.11.0</li>
</ul>

Note: All anaconda packages are listed in the provided environment.yml for ease of installation.

Be sure to edit the environment.ylm file and update the prefix directory at the very end of the file to be the absolute path where your anaconda3 environments are stored followed by the desired directory name for the new enviornment at the end.

This updated environment can be easily loaded by:
```
conda env create -f environment.yml
```
This environment.yml will install all necessary packages. If this step is taken, please skip the steps below.
However, if manual installation is required or you wish to use a different linear solver asides from the default MUMPS solver, please follow the steps below.

Please note that the manual steps are only validated for Linux systems. Installation of certain packages via Windows system during steps below may result in different package requirements as well as copying ipopt.exe to another location. If any issues are encountered please refer to IPOPT and Pyomo documentation directly for details and troubleshooting.

# Installation steps

$ENV is the python environment directory of your choice. This is where you create your python environment. For instructions below, replace the string $ENV with this directory using the absolve path
<ol>
<li>Create conda environment</li>
    <ul><li>&lt;conda create –prefix=$ENV python=3.8.5&gt;</li></ul>
<li>Activate conda environment</li>
    <ul><li>&lt;conda activate $ENV&gt;</li></ul>
<li>Install pyomo 
    <ul><li>&lt;conda install -c conda-forge pyomo==6.4.0&gt;</li></ul>
<li>Install pandas
    <ul><li>&lt;conda install pandas==1.2.5&gt;</li></ul>
<li>Install openpyxl
    <ul><li>&lt;conda install openpyxl==3.0.7&gt;</li></ul>
<li>Install cobra (conda does not install cobra due to python version not meeting requirements)
    <ul><li>&lt;pip install cobra==0.26.3&gt;</li></ul>
<li>Install scipy package
    <ul><li>&lt;conda install -c conda-forge scipy==1.8.0&gt;</li></ul>
<li>Install IPOPT 3.14.0
    <ul><li>&lt;conda install -c conda-forge ipopt==3.14.0&gt;</li></ul>
</ol>

# Optional but recommended: Harwell Subroutine Libraries
The above protocol will install pyomo and related dependencies to run the toy kinetic model script. However, the automatic installation of ipopt through conda-forge will install the MUMPS linear solver. If you wish you use a different solver, such as the Harwell Subroutine Libraries (HSL), instead of executing step 8, do the following:
<ol>
<li>Apply for an academic license for HSL package at https://www.hsl.rl.ac.uk/</li>
    <ul><li>Download file {coinhsl-2021.05.05.tar.gz}</li></ul>
<li>Download IPOPT 3.14.0 at https://github.com/coin-or/Ipopt/releases/tag/releases%2F3.14.0</li>
    <ul><li>Download source code {Ipopt-releases-3.14.0.tar.gz}</li></ul>
<li>Change to $ENV directory</li>
    <ul><li>&lt;cd $ENV&gt;</li></ul></li>
<li>Clone and install ThirdParty-HSL, this package helps build and install the HSL routines for IPOPT</li>
    <ol>
    <li>&lt;git clone https://github.com/coin-or-tools/ThirdParty-HSL.git -b stable/2.1&gt;</li>
    <li>&lt;cd ThirdParty-HSL&gt;</li>
    <li>Move file {coinhsl-2021.05.05.tar.gz} into this folder</li>
    <li>Unpack the HSL files</li>
        <ol><li>&lt;gunzip coinhsl-2021.05.05.tar.gz&gt;</li>
            <li>&lt;tar xf coinhsl-2021.05.05.tar&gt;</li></ol>
    <li>Rename directory and move to directory</li>
        <ol><li>&lt;ln -s coinhsl-2021.05.05 coinhsl&gt;</li>
            <li>&lt;cd coinhsl&gt;</li></ol>
    <li>Ensure lapack is loaded</li>
        <ul><li>&lt;module load lapack&gt;</li></ul>
    <li>Set configure file</li>
        <ul><li>&lt;./configure -prefix=$ENV&gt;</li></ul>
    <li>Build libraries, install libraries, and add environment variable</li>
        <ol><li>&lt;make&gt;</li>
            <li>&lt;make install&gt;</li>
            <li>&lt;libtool finish $ENV&gt;</li></ol>
    <li>Configure in ThirdParty-HSL as well</li>
        <ol><li>&lt;cd ThirdParty-HSL&gt;</li>
            <li>&lt;./configure -prefix=$ENV&gt;</li>
            <li>&lt;make&gt;</li>
            <li>&lt;make install&gt;</li></ol>
   </ol>
<li>Clone and install ThirdParty-ASL, this package creates an executable ipopt when compiling. This executable is required for pyomo to run IPOPT.</li>
    <ol><li>Go back to $ENV directory</li>
        <ul><li>&lt;cd $ENV&gt;</li></ul>
    <li>&lt;git clone https://github.com/coin-or-tools/ThirdParty-ASL -b stable/2.0&gt;</li>
    <li>&lt;cd ThirdParty-ASL&gt;</li>
    <li>&lt;./get.ASL&gt;</li>
    <li>&lt;./configure -prefix=$ENV&gt;</li>
    <li>&lt;make&gt;</li>
    <li>&lt;make install&gt;</li>
    </ol>
<li>Compile and install IPOPT</li>
    <ol>
    <li>Go back to $Env directory</li>
        <ul><li>&lt;cd $ENV&gt;</li></ul>
    <li>Move file {Ipopt-releases-3.14.0.tar.gz} into current directory</li>
    <li>Unpack files</li>
        <ol><li>lt;gunzip {Ipopt-releases-3.14.0.tar.gz}</li>
            <li>&lt;tar xvf {Ipopt-releases-3.14.0.tar}</li></ol>
    <li>Rename directory</li>
        <ul><li>&lt;mv {Ipopt-releases-3.14.0} Ipopt&gt;</li></ul>
    <li>Change to Ipopt directory</li>
        <ul><li>&lt;cd Ipopt&gt;</li></ul>
    <li>Create build directory and move to build directory</li>
        <ol><li>&lt;mkdir build&gt;</li>
            <li>&lt;cd build&gt;</li></ol>
    <li>Configure the build with openmp flags, the OpenMP flag allows for MA86 run multiple-threads</li>
        <ul><li>(if compiler is mkl may need to module load mkl first)</li>
            <ul>&lt;../configure --prefix=$ENV ADD_CFLAGS=-fopenmp ADD_FFLAGS=-fopenmp ADD_CXXFLAGS=-fopenmp&gt;</li></ul>
        </ul>
        <li>Make,test, and build</li>
        <ol><li>&lt;make&gt;</li>
            <li>&lt;make test&gt;</li>
            <ul><li>With HSL solver in place, all test should pass</li></ul>
            <li>&lt;make install&gt;</li></ol>
    </ol>
<li>Add PATH to the ipopt executable</li>
        <ul><li>export PATH=$PATH:$ENV/Ipopt/build/src/Apps/AmplSolver/</li></ul>

</ol>

# Optional: PARDISO Solver
PARDISO is another solver that could also be used within Pyomo with KETCHUP as part of IPOPT. Instructions for this optional solver are the following.
<ol>
<li>Download pardiso license and library package from https://pardiso-project.org/</li>
    <ol>
    <li>Download libpardiso600_GNU720-X86-64.so</li>
    <li>You will be given a license key in the email when applying. Save the license key as {pardiso.lic} and copy the key inside</li>
    <li>Place the lic file in the directory</li>
    </ol>
<li>Compile and install IPOPT</li>
    <ol>
    <li>Go back to $Env directory</li>
        <ul><li>&lt;cd $ENV&gt;</li></ul>
    <li>Move file {Ipopt-releases-3.14.0.tar.gz} into current directory</li>
    <li>Unpack files</li>
        <ul><li>&lt;gunzip {Ipopt-releases-3.14.0.tar.gz}</li>
            <li>&lt;tar xvf {Ipopt-releases-3.14.0.tar</li></ul>
    <li>Rename directory</li>
        <ul><li>&lt;mv {Ipopt-releases-3.14.0} Ipopt&gt;</li></ul>
    <li>Change to Ipopt directory</li>
        <ul><li>&lt;cd Ipopt&gt;</li></ul>
    <li>Create build directory and move to build directory</li>
        <ul><li>&lt;mkdir build&gt;</li>
            <li>&lt;cd build&gt;</li></ul>
    <li>Configure the build with openmp flags, the OpenMP flag allows for MA86 run multiple-threads</li>
        <ul><li>&lt; ../configure --prefix=$MYENV ADD_CFLAGS=-fopenmp ADD_FFLAGS=-fopenmp ADD_CXXFLAGS=-fopenmp --with-pardiso="/gpfs/group/cdm8/default/jack/pyomo/4/pyomo_4/libpardiso600-GNU720-X86-64.so" &gt;
            <li>ADD_CXXFLAGS=-fopenmp&gt;</li></ul>
    <li>Make,test, and build</li>
        <ul><li>&lt;make&gt;</li>
            <li>&lt;make test&gt;</li>
                <ul><li>Tests will not pass, few linking to do after</li></ul>
            <li>&lt;make install&gt;</li></ul>
    </ol>
<li>Link LD_PRELOAD variable from where libblas and liblapack is installed.</li>
    <ol><li>The information should be located in {libipopt.so}
        <li>Change directory to find {libipopt.so} usually stored at:
            <ul><li>&lt;cd $ENV/Ipopt/build/src/.libs/&gt;</li></ul>
        <li>Look at installation directories
            <ul><li>&lt;ldd libipopt.so&gt;</li></ul>
        <li>Set LD_PRELOAD path to the {libblas.so.3} and {liblaplack.so.3} files
            <ul><li>&lt;export LD_PRELOAD=”/lib64/libblas.so.3:/lib64/liblapack.so.3”&gt;</li></ul>
    </ol>
</ol>

# Additional notes
If you have already ran &lt;make&gt;, &lt;make install&gt; and wish to restart with a clean directory, use the following commands in that directory:
```
make uninstall
make distclean
```

To run with multiple threads, set the environment variable OMP_NUM_THREADS before running python
For example, to run with 24 threads use the command:
```
OMP_NUM_THREADS=24
```

To properly run the ipopt solver under ma86 you need to create a file {ipopt.opt} and write in the file &lt;linear_solver ma86&gt;.
For PARDISO, use
```
linear_solver pardiso
```

Pynumero installation use
```
pyomo download-extensions
pyomo build-extensions
```


