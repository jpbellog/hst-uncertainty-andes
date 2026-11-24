# Inputs

Starting files of the reconstruction. Their provenance (sources, sha256 prefixes and readers) is set
out in `rotation_models/README.md`, section "Inputs".

| Folder | Content |
|---|---|
| `rotations/Global_Model_R_BG_2025.rot` | Starting rotation model; the source of each row is in its comment |
| `deposits/` | Porphyry (`Porphyry_Plate_ID_geodyn_v7.csv`) and epithermal (`Epithermal_Plate_Id_geodyn_v6_thesis.csv`) deposit databases of this work |
| `geochemistry/Full_Data_v11_tt_with_Clusters_v12_geodyn2.xlsx` | Dated igneous samples (magmatic background) of this work |
| `margin/` | Restored trench (`Margin_SAM_0.1.csv`, trench nodes every 0.1 Myr) and the present-day trench (`Margin_SAM.shp` and companions) |
| `validation/` | Files of the validations of Supplementary Section S5: GPlates motion paths (`motion_paths/`, 10, 36 and 56 Ma, f_PAC 0.36, 0.44 and 0.5) and reconstructed tracks (`HST_Margin/`, 0 to 70 Ma every 5 Myr) |

The deposit and sample databases update the compilation published with Bello-González et al. (2026),
doi:10.1016/j.pce.2025.104202. The scripts read these files and never write to this folder.
