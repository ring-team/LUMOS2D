# LUMOS2D
LUMOS regroups python algorithms to locally update 2D triangular meshes. 

### Requirement ###

* at least the CMake-3.0.0
* GCC 4.8 or newer
* **Compiled** MMG (https://github.com/MmgTools/mmg)
* **Compiled** VTK (https://www.vtk.org/files/release/9.0/VTK-9.0.1.tar.gz)
* python version 3.6 or 3.7
* gmsh library

### Configuring LUMOS ###

#### Linux ####
Obtain VTK (only to run the jupyter notebook):
Download VTK : https://www.vtk.org/files/release/9.0/VTK-9.0.1.tar.gz

Install VTK
```
tar -xvzf VTK-9.0.1.tar.gz
cd VTK-9.0.1 && mkdir build && cd build
cmake ..
make && sudo make install 
```

If OpenGL library is not found:
```
sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
```

Obtain MMG (for more information : https://github.com/MmgTools/Mmg/wiki/Setup-guide):

```
git clone https://github.com/MmgTools/mmg
cd mmg
git checkout develop
```
```
mkdir build
cd build && cmake ..
```

If VTK is correctly detected, in the terminal this line should appear:
```
-- Compilation with VTK: add vtp and vtu I/O.
```

If VTK is not correctly detect you can test:
``` 
cmake -D VTK_DIR=<path-to-vtk-install-dir> ..
```

Compilation:
```
make [-j]
```

Obtain gmsh and shapely
```
pip install gmsh-api
pip install Shapely (for signed distance computation)
```

Obtain jupyter:
```
pip install jupyterlab
```

Obtain LUMOS2D
```
git clone https://github.com/ring-team/LUMOS.git
```

#### Windows ####

As branch develop of Mmg code is not always steady on Windows, we chose to keep the code for linux use only.


### Launch Jupyter-notebook ###
Install some dependencies for visualization:
```
pip install scipy
pip install matplotlib
pip install meshio
pip install vtk
```

```
cd LUMOS/python
jupyter-notebook
```
Double click on paper_results.ipynb

Change the mmg_path at the end of the first cell with your path to mmg2d_O3 executable file.

To run a cell, press Shift+Enter 
