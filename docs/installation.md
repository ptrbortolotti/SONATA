# Installation
How to install OpenMPI, mpi4py, PETSc and petsc4py for use with OpenMDAO and SONATA on Ubuntu with Conda.

### 1. Prerequisites
  - Ubuntu:
  
    ```
    sudo apt-get install libibnetdisc-dev
    sudo apt-get install libblas-dev libatlas-dev liblapack-dev
    ```

### 2. Install OpenMPI 1.8.3

  - Download OpenMPI: https://www.open-mpi.org/software/ompi/v1.8.3/
  
  - Extract and configure
      ```
      cd ~/Downloads/openmpi-1.8.3
      mkdir build
      cd build
      # I prefer to keep /usr/local clean, so I put it in /opt/openmpi
      ./configure --prefix=/opt/openmpi --with-devel-headers --enable-binaries
      make
      sudo make install
      ```
      
  - Add the following to your bash profile  
      ```
      export LD_LIBRARY_PATH=LD_LIBRARY_PATH=/opt/openmpi/lib:$LD_LIBRARY_PATH
      export PATH=/opt/openmpi/bin:$PATH
      ```


### 3. Install an up-to-date anaconda python=3.7 distribution

### 4. Create a sonata environment with the full basic anaconda modules and activate it:
    ```
    $ conda create -n sonata anaconda
    $ source activate sonata
    ```

### 5. Install mpi4py
    Again, make sure you have a valid openmpi (1.8.3 in my case) version
    installed and the $PATH and $LD_LIBRARY_PATH Set correctly:

    ```
    git clone https://github.com/mpi4py/mpi4py.git ./mpi4py.git
    cd mpi4py.git
    python setup.py build --mpicc=/opt/openmpi/bin/mpicc               
    pip install .
    ```

### 6. Install PETSc and petsc4py
    ```
    git clone https://github.com/petsc/petsc.git ./petsc.git
    cd petsc.git
    
    ./configure
    python setup.py build
    python setup.py install --record files.txt
    cd ..
    git clone https://bitbucket.org/petsc/petsc4py.git ./petsc4py.git
    cd petsc4py.git
    python setup.py build
    python setup.py install --record files.txt
    ```

### 4. Install the openmdao module
    ```
    $ conda install -c conda-forge openmdao
    ```

### 5. Install the pythonocc-package
    ```
    $ conda install -c conda-forge pythonocc-core=7.4.0
    ```
    
### 6. Install the pint module
	This is used to change units in the SONATA/CBM - DYMORE interface.
    ```
    conda install -c conda-forge pint
    ```

### 7. Install the intervaltree package
	This is used for structuring the topology and the calculation of layup coordinates.
    ```
    $ conda install -c conda-forge intervaltree
    ```

### 8. Install the shapely package
	This is used for the discretization and approximation of offset curves during the topology generatio
    ```
    $ conda install -c conda-forge shapely
    ```

### 9. Install the triangle package
	This is used for the unstructured triangulation of the cavities and balance weight materials during the meshing process:
    ```
    $ pip install triangle
    ```

## References
For more instructions see:
 - https://gist.github.com/mrosemeier/088115b2e34f319b913a#file-mpi_and_petsc_for_openmdao_osx_ubuntu_conda-md
 
