# lightning_line_charges
MATLAB suite that solves for charge transfer in a lightning flash using LMA/LDAR and E-Field data. To find
the charge transfer, we create lines of best fit corresponding to regions of charge and their
images, identified by the LMA/LDAR data. The lines are produced using a dimension reduction technique
from Principal Components Analysis on the coordinates in 3 space from the LMA/LDAR data. Then, we integrate
over the line to produce a system of equations which is solved for each E-Field data source, yielding 
the charge transferred in the lightning flash. 

LCD_MAIN.m - The main script file.

intergrl.m - Function for integrating the best fit lines.

princom.m  - Function for performing principal components analysis on the LMA/LDAR data.

field.txt  - User entered delta - E values for each field mill.

mills.txt  - Coordinates for each field mill at Kennedy Space Center.
