# SONATA
<img src="docs/img/logo_wframe.png" align="left"  width="140">
'SONATA' is a toolbox for Multidiciplinary Rotor Blade Design Environment for Structural Optimization and Aeroelastic Analysis.
Structural helicopter rotor blade optimization comprises classical aeroelastic problems, where the aerodynamic behavior, the structural elasticity and vibrational dynamics have to be studied simultaneously. Since the dynamic and modal behavior is strongly related to the structural properties of the rotor blades, adjusting these properties is essential for an effective optimization. Nevertheless, identifying constraints based on elemental matrices to keep the solution within feasible boundaries is often a protracted and iterative task. The herein presented definition of the rotor blade topology is deliberately associated to the production of composite rotor blades. Thus, manufacturability is inherent from the geometric layup definition. Using orthogonal projection with corner-style differentiation the cross-section is discretized and processed by the Variational Asymptotic Beam Sectional Analysis (VABS) afterwards. [more](docs/intro.md)

## Installation

To use the full functionality of SONATA a bunch of installations have to be made and packages to be gathered. In this section a brief installation guide is presented that will help the user to install it properly. 
SONATA is developed to work with a python version >3.6. An old python 2.7 release can be found under the tag v0.1


1. A **python > 3.6** distribution is needed. It is recommended to use anaconda for easier package management https://www.anaconda.com/download/

2. Install the **pythonocc** precompiled binaries for MacOSX/Linux/Windows 32 or 64 with the amazing conda package management system. Simply run the following commands in the terminal.
   ```	conda install -c conda-forge -c dlr-sc -c pythonocc -c oce pythonocc-core==0.17.3	```

3. Install the **pint** module. This is used to change units in the SONATA/CBM - DYMORE interface.
   ``` conda install -c conda-forge pint ```

4. Install the **intervaltree** package. This is used for structuring the topology and the calculation of layup coordinates. 
   ``` conda install -c conda-forge intervaltree ```

5. Install the **shapely** package. This is used for the discretization and approximation of offset curves during the topology generation process: *Windows*: Install the precompiled binaries from the /package directory by running the following command: ```pip install Shapely-1.6.4.post1-cp36-cp36m-win_amd64```; *Linux*: ```pip install shapely==1.6.4```

6. Install the **triangle** package. This is used for the unstructured triangulation of the cavities and balance weight materials during the meshing process: *Windows*: Install the precompiled binaries from the /packages directory by running the following command: ```pip install packages/triangle-20170106-cp27-cp27m-win_amd64.whl```; *Linux*: ```pip install triangle```

7. Install the **openmdao** package. This is the python package that provides the necessary framework for SONATA. you can either use the pip to install the openmdao or clone it directly from https://github.com/OpenMDAO/OpenMDAO ```pip install openmdao```

8. Test the installation and all packages by excecuting the folloging python script:
   ```	python test_install.py```	

9. Now you can download or clone the repository and execute the main SONATA script. 
   ```	python SONATA.py```

## Usage



## Developers - Guide

please read the [developers-guide](docs/developer-guide.md)

* [PythonOCC](http://www.pythonocc.org/)
* [openMDAO](http://openmdao.org/)

* [OpenCascadeTechnology Documentation](https://www.opencascade.com/doc/occt-6.9.1/refman/html/index.html)
* [PythonOCC API Documentation](http://api.pythonocc.org/)
* [OpenMDAO Documentation](http://openmdao.org/twodocs/versions/latest/)
* [Dymoresolutions User's Manual](http://www.dymoresolutions.com/UsersManual/UsersManual.html)
* [Numpydoc docstring guide](https://numpydoc.readthedocs.io/)


## Publications:

**Pflumm, T., Garre, W., Hajek, M.:** A Preprocessor for Parametric Composite Rotor Blade Cross-Sections, 44th European Rotorcraft Forum, Delft, The Netherlands, 2018  [[pdf]](docs/Pflumm,%20T.%20-%20A%20Preprocessor%20for%20Parametric%20Composite%20Rotor%20Blade%20Cross-Sections%20(2018,%20ERF).pdf) [[more…\]](https://mediatum.ub.tum.de/604993?query=Pflumm&show_id=1455385) [[BibTeX\]](https://mediatum.ub.tum.de/export/1455385/bibtex)

**Pflumm, T., Rex, W., Hajek, M.:** Propagation of Material and Manufacturing Uncertainties in Composite Helicopter Rotor Blades, 45th European Rotorcraft Forum, Warsaw, Poland, 2019 [[more…\]](https://mediatum.ub.tum.de/1520025) [BibTeX\]](https://mediatum.ub.tum.de/export/1520025/bibtex)


## Acknowledgment:

This work was supported by the German Federal Minisfor Economic Affairs and Energy through the German Aation Research Program LuFo V-2 and the Austrian Rsearch Promotion Agency through the Austrian Research Program TAKE OFF in the project VARI-SPEED.

<img src="docs/img/acknowledgment.png" width="400">
