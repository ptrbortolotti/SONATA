

<img src="docs/logo_wframe.png" align="left"  width="140">

# SONATA: 

Multidiciplinary Rotor Blade Design Environment for Structural Optimization and Aeroelastic Analysis

## **Abstract:**

Structural helicopter rotor blade optimization comprises classical aeroelastic problems, where the aerodynamic behavior, the structural elasticity and vibrational dynamics have to be studied simultaneously. Since the dynamic and modal behavior is strongly related to the structural properties of the rotor blades, adjusting these properties is essential for an effective optimization. Nevertheless, identifying constraints based on elemental matrices to keep the solution within feasible boundaries is often a protracted and iterative task. The herein presented definition of the rotor blade topology is deliberately associated to the production of composite rotor blades. Thus, manufacturability is inherent from the geometric layup definition. Using orthogonal projection with corner-style differentiation the cross-section is discretized and processed by the Variational Asymptotic Beam Sectional Analysis (VABS) afterwards.



[TOC]

## Introduction:

The large number of constraints and design drivers from various disciplines makes the helicopter rotor blade development process difficult, time consuming and costly.
The entire design process represents a classical aeroelastic problem, where the aerodynamic behavior, the structural elasticity and vibrational dynamics have to be studied simultaneously. The behavior can therefore not be examined with separate analysis of the different disciplines \cite{Tarzanin1998}. The integration of all the appropriate disciplines in the design process implies not only limitations on the design from various disciplines, but also defining and accounting for interactions so that the disciplines influence design decisions simultaneously rather than sequentially \cite{Adelman1989}.
Historically, the design and development of improved or entirely new rotor blades is conducted by departments in a company that maintain their separate simulation codes for performing their specific tasks. The aerodynamics department is responsible for performance calculations, aero-acoustics, rotor-wake interaction, unsteady airload prediction and computational fluid dynamics while the dynamics department focuses on rotor vibratory loads, stability and aeroelastic models\cite{Tarzanin1998}. The structural department determines the elastic properties as well as strength and fatigue characteristics.  A Blade and Rotor Design Department often bundles the different aspects while considering materials, manufacturability, maintainability and safety requirements. \cite{Tarzanin1998}
This modular approach narrows the scope of solutions, because each department focuses on individual objectives satisfied by individual design parameters. Mutual interactions can only be covered by numerous iterations.
In contrast to that, a multidisciplinary approach offers a more systematic development process that is able to design a better helicopter rotor \cite{Adelman1989}. Because of the impact the rotor behavior has on the overall performance of the helicopter and on customer noticeable vibratory characteristics, rotor aeroelastic effects should be considered in the earliest stages of the design process \cite{Rohl2012a}.
In the last 25 years, researchers have repeatedly stated the need for a design methodology and optimization framework that combines computational efficiency of a beam description in aeromechanic analysis with a rotor blade structural model that is capable at describing realistic composite rotor blade cross-sections with respect to the structural properties, applied load, stress and strain distributions as well as design constraints \cite{Friedmann1991, Weller1988, Lim2016}. 

## Framework:

Our multidisciplinary rotor blade design framework is named SONATA (**S**tructural **O**ptimization a**n**d **A**eroelas**t**ic **A**nalysis) and is illustrated in the following figure \ref{fig:framework}. Like most environments it comprises of **three** main components that are wrapped into an optimization framework. 

![test123](docs/environment.png)

Fig. 1: SONATA: Multidisciplinary Rotor Blade Design Environment for Structural Optimization and Aeroelastic Analysis embedded in OpenMDAO.

1. As a **first** component, the current state of the art involves an aeromechanical analysis of rotorcraft blades which includes flexible multibody dynamics, nonlinear finite elements and various rotorcraft aerodynamic models. They are often referred to as Comprehensive Analysis. Examples are the widely used Comprehensive Analytical Model of Rotorcraft Aerodynamics and Dynamics II (CAMRAD II) \cite{Johnson2013} and the software Dymore \cite{Bauchau2001} beyond several others. Both of these codes are presently in use in the rotorcraft industry, academic institutions and government laboratories. The quality of the predictions have been documented in numerous publications. In our SONATA environment Dymore was chosen as aeromechanic tool for both a dynamic analysis in the time domain as well a modal analysis within the frequency domain. In this context classical 1D-beam elements are used to describe the rotor blade due to the much simpler mathematical formulation and reduced computational effort compared to a full three-dimensional finite element model of the composite rotor blade \cite{Datta2011a}. Typically, this approach decouples the realistic composite blade definition and the manufacturability constraints from the aeromechanic analysis and the predesign of structural blade properties. That way, problems in the blade design cannot be discovered until later in the process where changes are costly and time consuming \cite{Rohl2012b}. 
2. Although the three-dimensional finite element method is the most accurate approach to model realistic rotor blades, it is still not appropriate for the use in rotor blade predesign \cite{Li2008, Datta2011a}. The slender characteristic of rotor blades allows the simplification to treat them as one-dimensional body \cite{Yeo2010}. Cesnik and Hodges \cite{Cesnik1995} formulated the Variational Asymptotic Beam Sectional Analysis (VABS) to accurately represent the behavior that is associated with the reduction of two-dimensions. In other words, this method splits the three-dimensional elastic problem into a two-dimensional linear cross-section analysis and a one-dimensional nonlinear beam analysis, which is able to consider initially twisted and curved, anisotropic, non-homogeneous materials to model general composite cross-sectional geometries \cite{Cesnik1995, Li2008}. VABS is the **second** component of our environment. In the last 20 years, VABS and its variations have become a popular tool in rotor blade predesign and multidisciplinary rotor design optimization and their accuracy and efficiency has been validated in numerous publications \cite{Cesnik1995, Cesnik2004, Rohl2012a}.
3. Consequently, most researches have developed individual parametric mesh generators for the cross-sectional analysis, that reduces their structural model to few design variables in the process. Such a preprocessor for parametric composite rotor blade cross-sections is referred to as *SONATA-CBM* in this framework. It is the **third** component of the SONATA environment. 

Last but not least, the tree components are managed by an environment where design variables and objectives can be defined, constraints to be applied and solvers to be launched. The **SONATA** framework uses [OpenMDAO](http://openmdao.org/)  \cite{Gray_2014, hwang_thesis_2015, Heath2013}, an open-source computing platform for system analysis and multidisciplinary optimization, written in Python. It allows the user to break down the structure of complex optimization tasks into a hierarchic manner while managing the numerical methods. A Python-based wrapper for Dymore has been developed to integrate the dynamic and modal analysis into the OpenMDAO-driven optimizations. Consequently \textit{SONATA-CBM} has been written in Python using the Python wrapper for the CAD-Kernel Opencascade [pythonOCC](http://www.pythonocc.org/ ). 

**Why Python?**

-  Python can be easy to pick up whether you're a first time programmer or you're experienced with other languages. 
-  Python is developed under an OSI-approved open source license, making it freely usable and distributable, even for commercial use. 
-  The Python Package Index (PyPI) hosts thousands of third-party modules for Python. Both Python's standard library and the community-contributed modules allow for endless possibilities. Two of the most important python modules used in SONATA are the openMDAO and the pythonocc module.
-  [openMDAO](http://openmdao.org/) is an open-source high-performance computing platform for systems analysis and multidisciplinary optimization, written in Python.
-   [pythonOCC](http://www.pythonocc.org/ )  is a python library whose purpose is to provide 3D modeling features. It is intended to developers who aim at developing CAD/PDM/PLM applications.

### 1. DYMORE (PYMORE):

this module comes [https://gitlab.lrz.de/wgarre/Pymore](https://gitlab.lrz.de/wgarre/Pymore)

### 2. VABS (Variational Asymptotic Beam Sectional Analysis):

Analyswift: Free Academic Liceses for Universities

### 3. SONATA-CBM:
*SONATA-CBM'*s composite topology generation originates from an arbitrary closed curve that can be obtained from various input formats that range from airfoil coordinate tables over a 3D CAD rotor blade surface  definition (.step or .iges) with radial station to a parameterized rotor blade with twist, planform, airfoil and chord-line distribution. In the case of the latter two, the 3D surface is intersected at a certain radial station to obtain once again a two-dimensional outer boundary of the cross-section. Figure \ref{fig:3dtopo} shows the resulting parameterized 3D surface of the UH-60A rotor blade with a cross-section topology at radial station $R=2000$mm.\\
While the following methodology is shown  with the example of the UH-60A rotor-blade, it should be noted that this procedure can be applied to any closed curve cross-section, and therefore be also used to model rotor blade root sections or any other composite beam cross-sections. 

<img src="docs/3dtopo.png" width="500">

<center>Fig. 2: Parameterized 3D surface of the UH-60A rotor blade created with twist, planform, airfoil and axis information from Davis [...] </center>

The process behind the composite topology generation is derived from the manufacturing process, where the layers are placed on top of each other in negative molds in a consecutive manner to avoid complex constraints in the optimization and to keep the solution within proper bounds. Each layer has an assigned material with start and end coordinates, a thickness and fiber orientation (see table \ref{tab:layup}). Every parameter or groups of them can serve as design variable in the later optimization. After the layup process on top of the outer boundary curve is completed, webs are introduced and subsequently new closed curved geometries are generated where the layup procedure is repeated. Cavities can be filled with core materials and additional trim masses can be inserted.
At first the outer boundary curve, represented as counterclockwise sets of consecutive B-splines, is defined in curve coordinates **s** between zero and one. The origin is typically located at the trailing edge (TE). The curve coordinate system propagates through the layers with an interval tree structure. It allows to efficiently find the intervals/layers that overlap and locate the corresponding coordinate for each layer. 
Subsequently, each layer is generated by the following consecutive steps. 

- Determine the relevant set of underlying B-Splines between \textit{start} and \textit{end} coordinate of the layer using an interval tree data structure.
- Discretize the set of B-Splines and perform an parallel offset to return an approximate representation of all points with a given thickness of each layer.
- Generate a new set of B-Splines by interpolation and add smooth layer cutoffs to connect the lower and upper set of B-Splines if necessary.

In table 1 the layup definition of the cross-section, illustrated in figure 2, is displayed. Note that the shown genetic composite cross-section of the UH-60A serves as demonstration of the modeling capabilities.



<img src="docs/layup.png" width="360">

<center>Table 1: Layup definition of figure 3 </center>

​						

## Installation

To use the full functionality of SONATA a bunch of installations have to be made and packages to be gathered. In this section a brief insallation guide is presented that will help the user to install it properly. 
SONATA is developed to work with a python version >3.6. An old python 2.7 release can be found under the tag v0.1

1. A python 2.7 distribution is needed. It is recommended to use use Anaconda for easier package management https://www.anaconda.com/download/
2. Install the **pythonocc** precompiled binaries for MacOSX/Linux/Windows 32 or 64 with the amazing conda package management system. Simply run the following commands in the terminal (for Windows users: execute the cmd command terminal):
    ```	conda install -c conda-forge -c dlr-sc -c pythonocc -c oce pythonocc-core==0.18	```

3. Install the **pint** module. This is used to change units in the SONATA/CBM - DYMORE interface.
    ``` conda install -c conda-forge pint ```

4. Install the **intervaltree** package. This is (will be) used for structuring the topology and the calculation of layup coordinates. 

  * ` conda install -c conda-forge intervaltree `

5. Install the **shapely** package. This is used for the discretization and approximation of offset curves during the topology generation process:
	* __Windows__: Install the precompiled binaries from the /package directory by running the following command: 
		
        ```pip install Shapely-1.6.4.post1-cp36-cp36m-win_amd64```
	* Linux: ```pip install shapely==1.6.4```
	
6. Install the **triangle** package. This is used for the unstructured triangulation of the core and balance weight materials during the meshing process:
	* __Windows__: Install the precompiled binaries from the /packages directory by running the following command: 
		
        ```pip install packages/triangle-20170106-cp27-cp27m-win_amd64.whl```
	* __Linux__: ```pip install triangle```

7. Install the **intervaltree** package. This is (will be) used for structuring the topology and the calculation of layup coordinates. 

  * ` conda install -c bioconda intervaltree `

8. Install the **openmdao** package. This is the python package that provides the necessary framework for SONATA. you can either use the pip to install the openmdao or clone it directly from https://github.com/OpenMDAO/OpenMDAO
	
	* `pip install openmdao`  
	* To use the pyoptsparse optimisation package within openmdao you need to install conda-	build ` conda install conda-build`. Then clone or download the repository from https://bitbucket.org/mdolab/pyoptsparse and build it like so ` conda build pyoptsparse `.   To use parallel computing features you need to follow the following instructions https://openmdao.readthedocs.io/en/1.7.3/getting-started/mpi_windows.html
	
9. Test the installation and all packages by excecuting the folloging python script:
	```	python test_install.py```

10. Now you can download or clone the repository and execute the main SONATA script. 
	```	python SONATA.py```


## Resources
* [PythonOCC](http://www.pythonocc.org/)
* [openMDAO](http://openmdao.org/)


#### Documentation for Developers:

* [OpenCascadeTechnology Documentation](https://www.opencascade.com/doc/occt-6.9.1/refman/html/index.html)
* [PythonOCC API Documentation](http://api.pythonocc.org/)
* [OpenMDAO Documentation](http://openmdao.org/twodocs/versions/latest/)

## Publications:
to follow...


## Referencencs:
[1] Tarzanin, F., and Young, D., “Boeing rotorcraft experience with rotor design and optimization,” 7th AIAA/USAF/ ASA/ISSMO Symp. Multidiscip. Anal. Optim., American Institute ofAeronautics and Astronautics, Reston,Virigina, 1998. doi:10.2514/6.1998-4733, URL [http://arc.aiaa.org/doi/abs/10.2514/6.1998-4733](http://arc.aiaa.org/doi/abs/10.2514/6.1998-4733).

[2] Adelman, H. M., and Mantay, W. R., “Integrated Multidisciplinary Optimization of Rotorcraft: A Plan for  Development,” Tech. rep., NASA, 1989.

[3] Rohl, P. J., Kumar, D., Dorman, P., Sutton, M., and Cesnik, C. E. S., “A Composite Rotor Blade Structural Design Environment for Aeromechanical Assessments in Conceptual and Preliminary Design,” American Helicopter Society 68th Annual Forum, American Helicopter Society, 2012. URL [http://ebooks.cambridge.org/ref/id/CBO9781107415324A009](http://ebooks.cambridge.org/ref/id/CBO9781107415324A009).

[4] Rohl, P., Dorman, P., Sutton, M., Kumar, D., and Cesnik, C., “A Multidisciplinary Design Environment for Composite Rotor Blades,” 53rd AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics and Materials Conference, American Institute of Aeronautics and Astronautics (AIAA), Reston, Virigina, 2012, pp. 1–15. doi:10.2514/6.2012-1842, URL [http://dx.doi.org/10.2514/6.2012-1842](http://dx.doi.org/10.2514/6.2012-1842). 

[5] Friedmann, P. P., “Helicopter Vibration Reduction Using Structural Optimization with Aeroelastic/multidisciplinary Constraints- A Survey,” Journal of Aircraft, Vol. 28, No. 1, 1991, pp. 8–21. doi:10.2514/3.45987, URL [http://dx.doi.org/10.2514/3.45987](http://dx.doi.org/10.2514/3.45987).

6] Weller, W. H., and Davis, M. W., “Wind Tunnel Tests of Helicopter Blade Designs Optimized for Minimum Vibration,”
American Helicopter Society 44th Annual Forum, 1988.

[7] Lim, J., Shin, S., and Kee, Y., “Optimization of Rotor Structural Design in Compound Rotorcraft with Lift Offset,” Journal of the American Helicopter Society, Vol. 61, No. 1, 2016, pp. 1–14. doi:10.4050/jahs.61.012005, URL [http://dx.doi.org/10.4050/JAHS.61.012005](http://dx.doi.org/10.4050/JAHS.61.012005).