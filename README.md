# GIPL2
====

Geophysical Institute Permafrost Laboratory 2 model version 01 (GIPL2v01) <br />
GIPL2 is a transient numerical model that employs phase changes and the effect of the unfrozen volumetric water content in the non-homogeneous soil texture. <br />
The original version of the model was developed by Romanovsky and Tipenko 2004 at the Geophysical Institute at the University of Alaska Fairbanks, and described in Marchenko et al., (2008). <br />
The present version has been significantly modified from its predecessor and using the IRF coding structure (see [Basic Machine Interface](http://csdms.colorado.edu/wiki/BMI_Description)). <br />
This version is maintained by Elchin Jafarov. <br />
Please cite Jafarov et al., (2012) when using the model. <br />

More details on model implementation could be found in Jafarov et al., (2012; 2014) and on the [CSDMS model site] (http://csdms.colorado.edu/wiki/Model:GIPL).

###References
Jafarov E.E., Nicolsky D.J., Romanovsky V.E., Walsh J.E., Panda S.K., Serreze M.C. 2014. The effect of snow: How to better model ground surface temperatures. Cold Regions Science and Technology, Volume 102, Pages 63-77, ISSN 0165-232X, [doi: 10.1016/j.coldregions.2014.02.007](http://www.sciencedirect.com/science/article/pii/S0165232X1400038X). <br />
Jafarov, E. E., Marchenko, S. S., and Romanovsky, V. E. 2012. Numerical modeling of permafrost dynamics in Alaska using a high spatial resolution dataset, The Cryosphere, 6, 613-624, [doi:10.5194/tc-6-613-2012](http://www.the-cryosphere.net/6/613/2012/tc-6-613-2012.pdf).

###Compile: 
**Windows**: Compile the gipl.f90 and call the executable file gipl.exe  <br />
**Linux**: Use Makefile to create executable, just type 'make' in your command line  <br />

###Run: 
Make sure to create a dump folder before running the executable file (cfg file in Setup below).  <br />

### Visualize:
The file with measured temperatures is 'mesres.txt'. The header for mesres is in the compare.m script. The command '>>compare(0)' executes the matlab script that plots the daily measured against calculated ground temperatures at four specified depth. 
![results](https://github.com/Elchin/GIPL/blob/master/results.png)

### Input/Output Setup:
Please see the 'cfg' file for more information on how to organize and input and output files. All input files should be stored in "in" folder. The path for the output files can be prescribed in the config file. Current output configuration includes 3 files: 'results.txt' file with daily output, 'mean.txt' with the yearly averaged data such as active layer thickness and freeze-up day, and the file 'start.txt' includes the temperature profile with depth for the last day of simulation. The result.txt has the following structure: time, upper_bnd_temperature, snow_depth, ground_temperatures, where 'ground_temperatures' assigned in grid.txt (see below). The mean.txt file has the same configuration and includes 3 more columns at the end. Everything in mean.txt is averaged yearly (see number_of_time_steps in config file below).

###Input data:
All input files are located in the "in" folder

**gipl_config.cfg** includes paths for input and uotput files and the correspoding setup paramteres. <br />
0/1: start from previous time step / start from the begining<br />
step | taum | tmin : <br />
    step is the timestep in the example it is 1<br />
    taum is the convergence parameter used by the stefan subroutine <br />
    tmin minimal timestep used in the Stefan subroutine <br />
begin | end : start and end, in the example it runs over one year from 0 to 1<br />
smoothing_factor | unfrozen_water_parameter | max number of iterations<br />
number_of_second_per_day [sec] | number_of_time_steps (in the example number of days in a year )<br />
sea_level | max_number_of_freezing_fronts [integer number]<br />
freezing_front_min | freezing_front_max depth [meters]<br />
saturation_coefficient (fraction of 1)<br />

File **input.txt** includes the total number of sites and the corresponding ids for the organic and mineral soils. In the current version, the number of sites is equal to 1.

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
vertical grid (in the example starts from the 1.5 meters above the ground up to 90 meters) 
minus sign corresponds to the values above the ground surface and plus corresponds to the values below the surface.
Let me denote the n elements of the grid as "gr"
The n+1 element of the grid corresponds to the number of output points (in the example it is 12)
The rest of the grid file correspond to indexes of the grid points (e.g. the number 40 below the 12 is the index of the 
gr(40)=0.001)

**initial.txt** includes all initial temperatures <br />
the code reads initial.txt file when in cmd.txt the first element is equal to 1. (e.i. start from the initial time step)
the first parameter in the file can be ignored / second parameter is the number of points (in the example it is 13)
The first column corresponds to the depth in meters and the second column to the temperature [Celsius] measured at that depth at time=0

**mineral.txt** is the most important file which includes all thermo-physical properties of the multilayered soil column
the first row can be ignored<br />
second row first element can be also ignored / second element in this row corresponds to the number of layer <br />
starting from row 3 to row 8 are thermo-physical properties of each layer.<br />
the first column is the volumetric water content (WVC) [fraction of 1]<br />
second and third columns are "a" and "b" coefficients of the unfrozen water curve (obtained from unfrozen water curve fitting) [dimensionless]<br />
forth and fifth columns are the thawed and frozen volumetric heat capacities [J/(m^3K)]<br />
six and seven columns are the thawed and frozen heat conductivities [W/(mK)]<br />
eighths column is the thickness of the corresponding layer<br />

File **organic.txt** has a similar structure to mineral and carries the parameters for the organic soil layer.
