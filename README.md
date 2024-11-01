
# UpsFrac
A modeling and upscaling software for two-dimensional fractured porous rocks.


>>> Perform DFM modeling using the ADFNE framework in the DFM_modeling folder.


>>> Subdivide the fractures into Cartesian grid blocks.


>>> Create input files for DFM simulations within Cartesian grid blocks.


>>> Conduct DFM simulations using the MRST framework in the DFM_simulation folder.


>>> Calculate the equivalent permeability (keq) for all Cartesian grid blocks.


>>> Visualize the equivalent permeability (keq).



UpsFrac
UpsFrac is an open-source software designed for modeling and upscaling permeability in fractured porous rocks, aimed at enhancing the quantitative and multi-scale analysis of fracture systems in geosciences.


Getting started
Download or clone this repository:
git clone git@github.com:chentao9330/UpsFrac.git

Once the repository has been downloaded or cloned, navigate your install of MATLAB to the checked out repository
cd UpsFrac;

Add the following folders to the path in MATLAB, and then you can start using the software

DFM_modeling folder based on MRST
DFM_simulation folder based on ADFNE


Software Usage Steps
step1： run   DiscreteFractureModeling.m , with the framework of ADFNE require: rndm_powerlaw.m; frac_determintic.txt
output: fracture geometric files, e.g., frac_0_1.txt

step2: run SubDivideToGrid.m
output: clipped fracture geometrie in coarse grid, e.g., resu_0_1.txt

step3: run DiscreteFractureInput.m, require: create InputTemplate_mpfa.m or use InputTemplate_tpfa.m
output: DFM simulation input files in folders, e.g., reali_1 folder

step4: run GridDFMSimu.m, do DFM simulation with the framework of MRST 
output: flow simulation results in folder \reali_1\RESULT

step5: calculate keq, run CalcuKeq_0902.m in folder \reali_1\RESULT,
output: upscaled permmeability in file, epcomf.txt

step6: visualization of keq, run plot_keq.m,plot_keq_hist.m,plot_keq_ellipse.m based on epcomf.txt
output: figures of equivalent permeability in coarse grid


# UpsFrac

UpsFrac is a scale-up software designed to enhance the analysis of fracture systems in geosciences. This tool provides advanced algorithms for processing and interpreting geochemical data from geothermal reservoirs.

## Why Use UpsFrac?

UpsFrac allows researchers and professionals in geosciences to:
- Improve the accuracy of geochemical analyses.
- Efficiently process large datasets.
- Gain insights into the behavior of fracture systems in geothermal environments.

- 

# UpsFrac
A modeling and upscaling software for two-dimensional fractured porous rocks.


>>> Perform DFM modeling using the ADFNE framework in the DFM_modeling folder.


>>> Subdivide the fractures into Cartesian grid blocks.


>>> Create input files for DFM simulations within Cartesian grid blocks.


>>> Conduct DFM simulations using the MRST framework in the DFM_simulation folder.


>>> Calculate the equivalent permeability (keq) for all Cartesian grid blocks.


>>> Visualize the equivalent permeability (keq).



# UpsFrac

UpsFrac is an open-source software designed for modeling and upscaling permeability in fractured porous rocks, aimed at enhancing the quantitative and multi-scale analysis of fracture systems in geosciences.

## Getting Started

To get started, download or clone this repository:

```bash
git clone git@github.com:chentao9330/UpsFrac.git
```

Once the repository has been downloaded or cloned, navigate to the repository in your MATLAB installation:

```bash
cd UpsFrac
```

Add the following folders to the MATLAB path, and then you can start using the software:

DFM_modeling: developed based on MRST
DFM_simulation: developed based on ADFNE


Software Usage Steps
Step 1: Run DiscreteFractureModeling.m
Requirements: rndm_powerlaw.m, frac_deterministic.txt
Output: Fracture geometric files (e.g., frac_0_1.txt)

Step 2: Run SubDivideToGrid.m
Output: Clipped fracture geometry in a coarse grid (e.g., resu_0_1.txt)

Step 3: Run DiscreteFractureInput.m
Requirements: Create InputTemplate_mpfa.m or use InputTemplate_tpfa.m
Output: DFM simulation input files in folders (e.g., reali_1 folder)

Step 4: Run GridDFMSimu.m
Perform DFM simulation with the framework of MRST
Output: Flow simulation results in the folder \reali_1\RESULT

Step 5: Calculate upscaled equivalent permeability, keq, for a coarse grid
Run CalcuKeq.m in the folder \reali_1\RESULT
Output: Upscaled permeability in the file epcomf.txt

Step 6: Visualization of upscaled equivalent permeability keq
Run plot_keq.m, plot_keq_hist.m, and plot_keq_ellipse.m based on epcomf.txt
Output: Figures of equivalent permeability in the coarse grid


License
This project is licensed under the GPL 3 License. See the LICENSE file for details.

Contact
For questions or feedback, please contact Tao Chen at chentao9330@gmail.com.
