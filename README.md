# UpsFrac
A modeling and upscaling software for two-dimensional fractured porous rocks.



step1： run   DiscreteFractureModeling.m , require: rndm_powerlaw.m;frac_determintic.txt   >>> Perform DFM modeling using the ADFNE framework in the DFM_modeling folder.
output:frac_0_1.txt

step2: run SubDivideToGrid.m   >>> Subdivide the fractures into Cartesian grid blocks.
output:resu_0_1.txt

step3: run DiscreteFractureInput.m, require: create InputTemplate_mpfa.m or use InputTemplate_tpfa.m   >>> Create input files for DFM simulations within Cartesian grid blocks.
output:reali_1 folder

step4: run GridDFMSimu.m   >>> Conduct DFM simulations using the MRST framework in the DFM_simulation folder.
output:reali_1\RESULT folder

step5: run CalcuKeq.m in folder \reali_1\RESULT   >>> Calculate the equivalent permeability (keq) for all Cartesian grid blocks.
output:epcomf.txt

step6: run plot_keq.m, plot_keq_hist.m, plot_keq_ellipse.m   >>> Visualize the equivalent permeability (keq).
