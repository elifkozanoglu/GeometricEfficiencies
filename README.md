# GeometricEfficiencies
Geometric efficiency for superNEMO OM-Bismuth source pairs

This project revolves around the calculation and visualization of geometrical efficiencies of main wall optical module-Bismuth source pairs in superNEMO experiment.
total_eff_visualization.cpp: creates 2D histogram visualization of (1) total geometrical efficiencies multiplied by activities of the sources as measured in July 1 2018 by Miro and (2) geometrical efficiencies associated with a singular source
plotting.cpp: generates a plot of optical modules (as boxes) and sources (as green dots) according to superNEMO geometry
function_visualization.cpp: generates a plot of the solid angle function for visualization
source_positions.txt.in: contains positions of the 42 Bismuth sources as provided by Filip's CalibrationModule
source_activity.txt.in: contains activities associated with each Bismuth source and their associated numbers as numbered by superNEMO (not explicitly used in the code)
