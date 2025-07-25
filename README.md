# __ML/MM setup__
This repository contains useful scripts to quickly set up QM/MM calculations within ChemShell with a level of theory of B3LYP/6-31G(d) for the QM part and AMBER (ff19SB) for the MM part. 

The setup uses the MD trajectory to extract precatalytic frames, according to a 'precatalytic definition', from which QM/MM scans are used to determine the energetic contribution to the free energy barrier, $\Delta G^‡$.

Besides QM/MM, scripts are also provided to set up ML/MM calculations using UMA as level of theory for the ML part and AMBER for the MM part.   
