# code-for-Size-resolved aerosol optics reveal ozone trade-offs for staple crop productivity
This repository contains all R scripts, data processing pipelines, and figure generation codes for the paper:
“Size-resolved aerosol optics reveal ozone trade-offs for staple crop productivity”
The study quantifies the impacts of fine-mode aerosols (fAOD) and ozone (O₃) on crop yield and calorie supply in China, and explores air pollution control strategies that co-optimize food security outcomes.

Code Structure and Naming Convention
The codebase follows a modular and numbered naming convention to ensure clarity and reproducibility:

00_ prefix: Initialization scripts, including package loading (00_loadPackages.R), model formulas (00_loadFormulas.R), and common functions (00_loadFunctions.R).
01_ prefix: Geographic masking and statistical data processing, such as crop area, production, and population statistics (01_prepare_masks_stats.R, etc.).
02_ prefix: Data cleaning, transformation, and integration of remote sensing products and air pollution metrics (02_preprocess_ozone_data.R, etc.).
03_ prefix: Model fitting and counterfactual simulations for crop response to pollutants (03_run_main_models.R, etc.).
04_ prefix: Generation of main and supplementary figures for publication, including maps, line plots, and surface plots (04_plot_fig2ab.R, 04_plot_fig4.R, etc.).
Each script is self-contained and includes appropriate comments to guide replication and adaptation.
All analyses were conducted under R version 4.2.3, ensuring compatibility and reproducibility across platforms.
