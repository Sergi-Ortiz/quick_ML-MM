# __ML/MM Workflow__

This document aims to explicitly showcase the usage of `mlmm-toolkit` to swiftly study enzymatic reactions (or homogeneous reactions in any other system) with ML/MM calculations. 

>This tutorial assumes familiarity with all the methods used and that `mlmm-toolkit` is either installed as a local conda environment or in a SLURM cluster with the same conda environment. 

> The installation guide of `mlmm-toolkit` and its installation in a SLURM machine can be found in [TODO](). 

---

## 0. System preparation

To perform a ML/MM study of any chemical reaction, we essentially need the same information as for QM/MM: the solvated ligand-enzyme complex, its structure, coordinates and topology. 

`quick-ML/MM` is designed to accelerate this workflow by
1. extracting precatalytic frames (according to a predefined structural metric) from a MD trajectory. 

2. trimming the solvent and parametrizing the resulting system.

3. copying the necessary files for each frame (structure, coordinates and topology). 

Given a specific frame, we will need three core components: the structure `complex.pdb`, its coordinates `complex.inpcrd` and its topology, `complex.prmtop`. 

> Currently, only AMBER parameter files are supported. CHARMM and GROMACS will be supported in the future. 


## 1. ML region definition

As in QM/MM, we need to define the ML region. It is important to note that the ML/MM embedding used by `mlmm-toolkit` is a **mechanical embedding**, and as such, MM polarization effects on the ML region are not accounted for. An electrostatic ML/MM embedding would perform better but foundational MLIP models currently in use, such as UMA, do not take charges as an input feature vector. 

> **Note.** The drawback of a mechanical embedding can be minimized trivially in the case of ML/MM by making the ML region larger, such that the atoms that would/could polarize the ML region are also included (obviously at the expense of computational performance). 

1. Define the ML region with the `def_ml_region` CLI helper. 
- `def_ml_region -r complex.pdb -c ligand.pdb -o ml_region_v2.pdb --radius_ml 3.5 --radius_het 3.5 --include_H2O false --exclude_backbone true -e -e 664 663 662 544 540 365 360 596 399`. The `-e` flag explicitly includes given residues. 

> Ensure that atom order, atom names, residue IDs and residue names are
identical to those in the full `complex.pdb`.

---


## 3. ML/MM workflow layout

Each frame is expected to have this structure

```text
frame
├── src
├── parm
│   ├── complex.pdb
│   ├── complex.inpcrd
│   ├── complex.prmtop
│   └── ml_region.pdb
├── job.yaml # job config file
└── job.sh   # bash script to run job
```

- `src` directory will contain important structures, such as optimized structures, scans, GSM strings, TS structures etc. 
- `parm` contains the essential reference files for a ML/MM simulation
- `job.yaml` corresponds to the configuration file for a particular type of job (optimization, scan, GSM, TS, IRC, etc.)
- `job.sh` corresponding SLURM runner to run a particular `job.yaml` job. 

The canonical pipeline consists of 10 `job.yaml` files. In the LOX example, `2_scan.yaml` is not performed (explained later why). 

| Stage | YAML file(s)              | Purpose                                             | 
|------:|---------------------------|-----------------------------------------------------|----------------------|
| (i)   | `1_opt.yaml`              | Initial optimization of the snapshot                | 
| (ii)  | `2_scan.yaml`             | 1‑D bond‑length scan → seed Reactant/Product        | 
| (iii) | `3_opt.yaml`, `4_opt.yaml`| Relax *reactant* and *product* obtained by scan     | 
| (iv)  | `5_gs.yaml`               | Minimum Energy Path search by Growing String Method | 
| (v)   | `6_tsopt.yaml`            | **PHG‑Dimer** transition state refinement           | 
| (vi)  | `7_irc.yaml`              | IRC propagation                                     | 
| (vii) | `8_opt.yaml`, `9_opt.yaml`| Endpoint relaxation → Final Reactant/Product        |
| (viii)| `10_energy_summary.yaml`  | Table and Figures of $\Delta E,G$                   | 


Each `.yaml` job can be initiated with 
```bash
mlmm          1_opt.yaml          # (i)
bond_scan     2_scan.yaml         # (ii)
mlmm          3_opt.yaml          # (iii‑R)
mlmm          4_opt.yaml          # (iii‑P)
mlmm          5_gs.yaml           # (iv)
ts_search     6_tsopt.yaml        # (v)
mlmm          7_irc.yaml          # (vi)
mlmm          8_opt.yaml          # (vii‑R)
mlmm          9_opt.yaml          # (vii‑P)
energy_summary 10_energy_summary.yaml   # (viii)
```

When running in a SLURM environment, all common folders (`src/`, `parm/` and the corresponding `job.yaml`) are copied to stratch, this maintains the same workflow as if it were a local conda environment. 

---

## 4. Core workflow
### 4.1.  optimization – `1_opt.yaml`

```bash
mlmm 1_opt.yaml
```

Perform an optimization of the precatalytic frame frozing atoms further than 15 Å from the OH cofactor.

> In preprocessing, each precatalytic frame is trimmed such that only solvent within 17 Å of the OH cofactor is included. The idea of frozing to 15 Å is such that there is a 2 Å buffer of frozen solvent between the mobile solvent and the void. 

The `frozen_atoms` list is obtained with `frozen_atoms.py` [todo](), which must be then copied to `1_opt.yaml`.

```bash
python frozen_atoms.py -i complex.pdb -o frozen_atoms.txt
```

The final structure is written to `frame/dump/1_opt/final_structure.xyz` and must be copied to `frame/src/1_opt_final_geometry.xyz` (automatic in the .sh runner). 

```bash
cp ./dump/1_opt/final_structure.xyz ./src/1_opt_final_geometry.xyz
```

The resulting optimized `.pdb` structure file can be obtained from the `.xyz` via 

```bash
xyz_geom2pdb -i ./src/1_opt_final_geometry.xyz -o ./src/1_opt_final_geometry.pdb -r ./parm/complex.pdb
```

### 4.2. Obtain a representative product structure

To obtain a product structure, depending on the reaction, one could perform a 1‑D scan along the forming bond. This can be done via 

```bash
bond_scan 2_scan.yaml
```
and visualizing the resulting trajectory with 

```bash
xyz_geom2pdb -i ./dump/2_scan/final_geometries.trj -o ./src/2_scan_path.pdb -r ./parm/complex.pdb
```

from which two representative structures are picked: the *reactant* (`2_initial_reac.pdb`) and *product* (`2_initial_prod.pdb`) nodes for the GSM-DE calculation. 

Note however that for some reactions a 1-D scan will not produce the desired product. In this case, the product structure should obtained manually with a structure editor. This product structure will be used in the GSM-DE calculation but will not be the final product structure (which will come from the IRC scan).  

For the ALOX15/5S-HETE example, the product can be obtained by deleting the CH bond, creating the OH bond, adjusting the OH bond distance and HOH angle. This is automatized with `product.py` which can be used via

```bash
python product.py -i ./src/1_opt_final_geometry.pdb -o ./src/2_initial_prod.pdb
```


### 4.3  Reactant / product optimizations – `3_opt.yaml`, `4_opt.yaml`

Perform optimizations for the *reactant* (`frame/src/2_initial_reac.pdb`) and *product* (`frame/src/2_initial_prod.pdb`) images. This step is the same as for **4.1.**. Resulting optimized structures should be visualized and checked for artifacts, as they are essential for the GSM calculations. 

```bash
# run the optimization of the reactant
mlmm 3_opt.yaml

# copy the optimized reactant structure to src
cp ./dump/3_opt_reac/final_structure.xyz ./src/3_opt_reac.xyz

# visualize the structure
xyz_geom2pdb -i ./src/3_opt_reac.xyz -o ./src/3_opt_reac.pdb -r ./parm/complex.pdb

# repeat for 4_opt.yaml
mlmm 4_opt.yaml
cp ./dump/4_opt_prod/final_structure.xyz ./src/4_opt_prod.xyz
xyz_geom2pdb -i ./src/4_opt_prod.xyz -o ./src/4_opt_rpdo.pdb -r ./parm/complex.pdb
```

### 4.4  GSM-DE calculation – `5_gsm.yaml`

An approximate TS structure can be obtained via a GSM-DE calculation

```bash
mlmm 5_gsm.yaml
```

Once the calculation has finished, the HEI and node energies / string energies can be extracted via

```bash
# obtain the energies as a csv and the HEI as an xyz
trj2fig -i ./dump/5_gsm/current_geometries.trj -o ./src/5_gsm.csv --output-peak ./src/5_gsm_HEI.xyz

# copy resulting .trj file to src
cp ./dump/5_gsm/current_geometries.trj ./src/5_gsm.trj

# obtain a .pdb for the HEI for visualization
xyz_geom2pdb -i ./src/5_gsm_HEI.xyz -o ./src/5_gsm_HEI.pdb -r ./parm/complex.pdb
```

The GSM should be visualized and evaluate if the HEI can be reasonably identified as a TS. 

### 4.5  Transition state search – `6_tsopt.yaml`

Starts from `5_gsm_HEI.xyz` and refines the TS with the
**PHG‑Dimer** algorithm.  

```bash
ts_search 6_tsopt.yaml
```

The optimized TS is saved as `./dump/dimer/final_geometry.xyz` and can be visualized via 

```bash
# copy the TS structure to src
cp ./dump/6_/final_structure.xyz ./src/3_opt_reac.xyz

# visualize the structure
xyz_geom2pdb -i ./src/3_opt_reac.xyz -o ./src/3_opt_reac.pdb -r ./parm/complex.pdb

# TODO
# see what to do with ./dump/6_ts_vib
```

Aditionally, a vibrational analysis of the TS is conducted and a single imaginary frequency should appear. The movement associated to this frequency should be visualized (should correspond to bond formation/dissociation motion). 


---
TODO CONTINUE ADDING THE TUTORIAL
### 4.6  IRC – `7_irc.yaml`

Propagates the TS downhill in both directions.  
The final frames (`./dump/irc/backward_last.xyz` → `7_irc_backward_last.xyz`, 
`./dump/irc/forward_last.xyz` → `7_irc_forward_last.xyz`) provide improved reactant/product guesses.

### 4.7  Endpoint optimizations – `8_opt.yaml`, `9_opt.yaml`

Relax the IRC endpoints to obtain the *final* reactant and product
geometries (stored as `./dump/opt4/final_geometry.xyz` and `./dump/opt5/final_geometry.xyz`).

### 4.8  Energy summary – `10_energy_summary.yaml`

Compute electronic + vibrational contributions for the final
reactant, TS and product:

```bash
energy_summary 10_energy_summary.yaml
```

This command prints $\Delta G$, $\Delta G^{\ddagger}$, $\Delta E$ and $\Delta E^{\ddagger}$ and writes figures of energy diagrams.

> All intermediate files live under `./dump/`.
> In this example, all input geometries are stored in `./coord/`.
---