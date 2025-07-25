# __ML/MM setup__
This repository contains useful scripts to quickly set up QM/MM calculations within ChemShell with a level of theory of B3LYP/6-31G(d) for the QM part and AMBER (ff19SB) for the MM part. 

The setup uses the MD trajectory to extract precatalytic frames, according to a 'precatalytic definition', from which QM/MM scans are used to determine the energetic contribution to the free energy barrier, $\Delta G^‡$.

Besides QM/MM, scripts are also provided to set up ML/MM calculations using UMA as level of theory for the ML part and AMBER for the MM part.   

## __Simulation structure__
QM/MM and ML/MM simulation have been set up with a specific file structure in mind. The setup process is automatic, and only the MD files have to be provided. An example simulation structure is provided in `sim_structure` directory, which includes the `1_frames` directory, where the MD data is located, and `src` directory, where the source scripts (for either QM/MM or ML/MM methods) must be included. 


## __QM/MM scripts__
Located in `qmmm` directory.



## __ML/MM scripts__
Located in `mlmm` directory.

