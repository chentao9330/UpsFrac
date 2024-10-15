# UpsFrac
A modeling and upscaling software for two-dimensional fractured porous rocks.


 >>> Perform DFM modeling using the ADFNE framework in the DFM_modeling folder.
> >> step1： run   DiscreteFractureModeling.m , require: rndm_powerlaw.m;frac_determintic.txt  
output:frac_0_1.txt

>>> Subdivide the fractures into Cartesian grid blocks.
>>> step2: run SubDivideToGrid.m   
output:resu_0_1.txt

>>> Create input files for DFM simulations within Cartesian grid blocks.
>>> step3: run DiscreteFractureInput.m, require: create InputTemplate_mpfa.m or use InputTemplate_tpfa.m   
output:reali_1 folder

>>> Conduct DFM simulations using the MRST framework in the DFM_simulation folder.
>>> step4: run GridDFMSimu.m   
output:reali_1\RESULT folder

>>> Calculate the equivalent permeability (keq) for all Cartesian grid blocks.
>>> step5: run CalcuKeq.m in folder \reali_1\RESULT   
output:epcomf.txt

>>> Visualize the equivalent permeability (keq).
>>> step6: run plot_keq.m, plot_keq_hist.m, plot_keq_ellipse.m   
