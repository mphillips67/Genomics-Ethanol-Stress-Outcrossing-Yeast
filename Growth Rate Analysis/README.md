# Growth Rate Analysis

This folder contains scripts and data to compare growth rate data for different treatments and ancestral population in control YPD media, 6% ethanol YPD, and 10% ethanol YPD.

Input files for these analyses can be found on the Dryad repository associated with this project but are also included in this repository for convenience. 

## Folder Contents

### Input files

This folder contains input growth rate data from experimental yeast populations:

- **Eth_48_hr_timepoint_0.txt** — Plate reader assay measuring growth in the ancestral population across all media types for 48 hours.  
- **Plain_ypd_2021.txt** — Plate reader assay measuring growth in the ancestral population and experimental populations during cycle 15 in plain YPD.  
- **Moderate_ypd_2021.txt** — Plate reader assay measuring growth in the ancestral population and experimental populations during cycle 15 in 6% ethanol YPD.  
- **High_ypd_2021.txt** — Plate reader assay measuring growth in the ancestral population and experimental populations during cycle 15 in 10% ethanol YPD.  

### Data Format

Columns in these files generally include:

1. **time** — Time (in minutes) at which OD600 absorbance was measured.  
2. **blank** — OD600 absorbance measurements for blank media wells.  
3. **experimental/ancestral populations** — OD600 absorbance measurements for each experimental and/or ancestral population at each timepoint.  

### Scripts

This folder contains scripts to:
- Analyze growth rate differences between treatments.
- Generate growth rate analysis figures for the associated manuscript including Figure 1, Supplementary Figure 1, and Supplementary Figure 3. 

