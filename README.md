# GIPL2
====

Geophysical Institute Permafrost Laboratory 2 model version 01 (GIPL2v01) <br />
GIPL2 is a numerical transient model that employs phase changes and the effect of the unfrozen volumetric  water content in the non-homogeniuos soil texture. <br />
Original version of the model developed by Romanovsky and Tipenko 2004 at the Geophysical Institute at Univeristy of Alaska Fairbanks, and described in Marchenko et al., (2008). <br />
Current version been significanlty modefied from its predicessor and using the IRF coding design (see [Basic Machine Interface](http://csdms.colorado.edu/wiki/BMI_Description)). <br />
This version is maintained by E. Jafarov at INSTAAR, CU Boulder. <br />
Please cite Jafarov et al., (2012) when using the model. <br />

More details on model implementation could be found in Jafarov et al., (2012; 2014) and on the [CSDMS model site] (http://csdms.colorado.edu/wiki/Model:GIPL).

###References
Jafarov E.E., Nicolsky D.J., Romanovsky V.E., Walsh J.E., Panda S.K., Serreze M.C. 2014. The effect of snow: How to better model ground surface temperatures. Cold Regions Science and Technology, Volume 102, Pages 63-77, ISSN 0165-232X, doi: 10.1016/j.coldregions.2014.02.007.  
Jafarov, E. E., Marchenko, S. S., and Romanovsky, V. E.: Numerical modeling of permafrost dynamics in Alaska using a high spatial resolution dataset, The Cryosphere, 6, 613-624, doi:10.5194/tc-6-613-2012, 2012

###What to do first: 
**Windows**: Compile the gipl.f90 and call the executable file gipl.exe  <br />
**Linux**: Use Makefile to create executable, just type 'make' in your command line  <br />
Please see the 'cfg' file for more information on how to orginize and input and output files. All input files should be stored in "in" folder. The path for the output files can be prescribed in the config file. Current output configuration include 3 files: 'results.txt' file with daily output, 'mean.txt' with the yearly averaged data such as active layer thickness and freezeup day, and the file 'start.txt' includes the temperature profilee with depth for the last day of a simulation. The file with measured temperatures is 'mesres.txt'. The compare(0) is a matlab script that polts the output daily file against measured data for four specified depths and estimates the mean average error between the measured and the calculated ground temperatures. 

###Input data:
All input files are located in the "in" folder

**gipl_config.cfg** includes paths for input and uotput files and the correspoding setup paramteres. <br />
0/1: start from previous time step / start from the begining<br />
step | taum | tmin : <br />
	step is the timestep in the example it is 1<br />
	taum is the convergence parameter used by the stefan subroutine <br />
	tmin minimal timestep used in the Stefan subroutine <br />
begin end : start and end, in the example it runs over one year from 0 to 1<br />
smoothing factor | unfrozen water parameter | max number of iterations<br />
number of second in a day [sec] | number of time steps (in the example number of days in a year )<br />
sea level | max number of freezing fronts [integer number]<br />
Freezing front min and max depth [meters]<br />
saturation coefficient (fraction of 1)<br />

File **input.txt** includes the total number of sites and the correspoding ids for the organic and mineral soils. In the current version the number of sites is equal to 1.

**bound.txt** include upper boundary condition (in the example it is an air temperature T [Celsius])<br />
First row is a number of points (in the example number of day)<br />
Column 1: Time parameter (day number)<br />
Column 2: Temperature (daily averaged temperature)

**snow.txt** include snow depth [meters] (in the example it is daily snow depth)<br />
First row is a number of points (in the example number of day)<br />
Column 1: Time parameter (day number)<br />
Column 2: Snow depth (daily averaged)

**rsnow.txt** include snow thermal conductivity [W/(mK)] (in the example it is daily snow conductivity)<br />
First row is a number of points (in the example number of day) <br />
Column 1: Time parameter (day number)<br />
Column 2: Snow conductivity (daily averaged)

**grid.txt** includes number of grid point (n), <br />
vertical grid (in the example starts from the 1.5 meter above the ground up to 90 meters) 
minus sign correspond to the values above the ground surface and plus corresponds to the values below the surface.
Let me denote the n elements of the grid as "gr"
The n+1 element of the grid correspond to the number of output points (in the example it is 12)
The rest of the grid file correspond to indexes of the grid points (e.g. the number 40 below the 12 is the index of the 
gr(40)=0.001)

**initial.txt** includes all initial temperatures <br />
the code reads initial.txt file when in cmd.txt the first element is equal to 1. (e.i. start from the initial time step)
first parameter in the file can be ignored / second parameter is the number of points (in the example it is 13)
First column corresponds to the depth in meters and the second column to the temperature [Celsius] measured at that depth at time=0

**mineral.txt** is the most important file which include all thermo-physical properties of the multilayered soil column
first row can be ignored<br />
second row first element can be also ignored / second element in this row correspond to the number of layer <br />
starting from row 3 to row 8 are thermo-physical properties of each layer.<br />
first column is the volumetric water content (WVC) [fraction of 1]<br />
second and third columns are "a" and "b" coefficients of the unfrozen water curve (obtained from unfrozen water curve fitting) [dimensionless]<br />
forth and fifth columns are the thawed and frozen volumetric heat capacities [J/(m^3K)]<br />
six and seven columns are the thawed and frozen heat conductivities [W/(mK)]<br />
eighths column is the thickness of the corresponding layer<br />

File **organic.txt** have a similar structure to mineral and carries the parameters for the organic soil layer
