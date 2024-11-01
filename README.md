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

DFM_modeling: developed based on ADFNE  
DFM_simulation: developed based on MRST


## Software Usage Steps
### Step 1: Run DiscreteFractureModeling.m
Requirements: DFM_modeling, frac_deterministic.txt (optional)
Input: fracture geometry parameters
Output: Fracture geometric files (e.g., frac_0_1.txt)  

### Step 2: Run SubDivideToGrid.m
Input: dimensions of a coarse grid
Output: Clipped fracture geometry in a coarse grid (e.g., resu_0_1.txt)  

### Step 3: Run DiscreteFractureInput.m
Requirements: input file template for dfm simulation on each coarse grid, InputTemplate_tpfa.m or InputTemplate_mpfa.m
Output: DFM simulation input files in folders (e.g., reali_1 folder)  

### Step 4: Run GridDFMSimu.m
Perform DFM simulation with the framework of MRST
Requirements: DFM_simulation 
Output: Flow simulation results in the folder \reali_1\RESULT  

### Step 5: Calculate upscaled equivalent permeability, keq, for a coarse grid
Run CalcuKeq.m in the folder \reali_1\RESULT  
Output: Upscaled permeability in the file epcomf.txt  

### Step 6: Visualization of upscaled equivalent permeability keq
Run plot_keq.m, plot_keq_hist.m, and plot_keq_ellipse.m based on epcomf.txt  
Output: Figures of equivalent permeability in the coarse grid  


## License
This project is licensed under the GPL 3 License. See the LICENSE file for details. The ADFNE and MRST used should adhere to their respective licensing terms.

## Contact
For questions or feedback, please contact Tao Chen at chentao9330@gmail.com.
