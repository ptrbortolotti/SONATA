# SONATA: A Preprocessor for Parametric Analysis and Design of Composite Beam Cross-Sections in a Multidisciplinary Rotor Design Environment

<img src="logo.png" align="left" hspace="20" vspace="6"> **SONATA** is a preprocessor for parametric analysis and design of composite beam cross-sections in a multidisciplinary rotor design environment. A helicopter rotor blade represents a classical aeroelastic problem, where the aerodynamic behavior, the structural elasticity and vibrational dynamics have to be studied simultaneously.  While a geometric definition of a rotorblade with CAD tools is simple, the transfer to a meshed cross-sectional representation may prohibit automated design optimization. Consequently, most researches have developed individual parametric mesh generators for the cross-sectional analysis, that reduces their structural model to few design variables in the process. SONATA represents such a preprocessor.
SONATA is written in python and is using for a lot of operations the Opencascade (CAD) kernel with its python wrapper (pythonocc). 


#SONATA helps the engineer to parameterize a closed composite rotor blade crossection with multiple spars. It is specifically designed to be suited for helicopter rotor blade crossections of the blade aerodynamic section and elastic blade root. SONATA combines visualization and 2D-Finite Element discretisation of the crossection. 

The first part of the software contains a parametric topology generator 
The topology is saved as a .pkl and can be reloaded
The second part generates a mesh upon the topology, the mesh can be exported into a VABS and SECTIONBUILDER conform PATRAN mesh file .ptr

More to come...
<img src="\img\bugless_meshing.png" hspace="20" vspace="6" width="600">

## Resources
* [PythonOCC](http://www.pythonocc.org/)

## Documentation for Developers:

* [OpenCascadeTechnology Documentation](https://www.opencascade.com/doc/occt-6.9.1/refman/html/index.html)
* [PythonOCC API Documentation](http://api.pythonocc.org/)
* [PythonOCC API Documentation](http://api.pythonocc.org/)




## Installing
1. A python 2.7 distribution is needed. It is recommended to use use Anaconda for easier package management https://www.anaconda.com/download/
2. You can install precompiled binaries for MacOSX/Linux/Windows 32 or 64 thanks to the amazing conda package management system. Simply run:

```html
conda install -c conda-forge -c dlr-sc -c pythonocc -c oce pythonocc-core==0.18
```

3. install the **shapely** package. This is used for the discretization and approximation of offset curves during the topology generation process:
.* Windows: Install the precompiled binaries from the /package directory by running the following command: 
	```html
	pip install Shapely-1.5.17-cp27-cp27m-win_amd64.whl
	``` 
.* Linux: TBD

3. install the **triangle** package. This is used for the unstructured triangulation of the core and balance weight materials during the meshing process:
.* Windows: Install the precompiled binaries from the /package directory by running the following command: 
	```html
	pip install triangle-20170106-cp27-cp27m-win_amd64.whl
	```
.* Linux: TBD

4. install the **openmdao** package. This is used for the unstructured triangulation of the core and balance weight materials during the meshing process:
.* Windows: Install the precompiled binaries from the /package directory by running the following command: 
	```html
	pip install openmdao
	```
.* Linux: TBD

5. install the **intervaltree** package.



```html
<script src="https://d3js.org/d3.v4.js"></script>
```


import OCC
import intervaltree 
import shapely
import openmdao