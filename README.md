GIPL
====

Geophysical Instutite Permafrost Lab (GIPL) model
Original version developed by Genadii Tipenko and Vladimir Romanovsky 
at the Geophysical Institute Permafrost Laboratory, University of Alaska Fairbanks (2004)
Modified by Elchin Jafarov (2012)

Input data:
All input files are located in "in" folder

bound.txt include upper boundary condition (in the example it is an air temperature T [Celsius])
First row is a number of points (in the example number of day)
Column 1: Time parameter (day number)
Column 2: Temperature (daily averaged temperature)

snow.txt include snow depth [m] (in the example it is daily snow depth)
First row is a number of points (in the example number of day)
Column 1: Time parameter (day number)
Column 2: Snow depth (daily averaged)

rsnow.txt include snow thermal conductivity [W/(mK)] (in the example it is daily snow conductivity)
First row is a number of points (in the example number of day)
Column 1: Time parameter (day number)
Column 2: Snow conductivity (daily averaged)

cmd.txt include 
0/1: start from previous time step / start from the begining
step | taum | tmin : 
     step is the timestep in the example it is 1
     taum is the convergence parameter used by the stefan subroutine 
     tmin minimal timestep used in the stefan subroutine 
begin end : start and end, in the example it runs over one year from 0 to 1
smoothing factor | unfrozen water parameter | max number of iterations
number of second in a day [sec] | number of time steps (in the example number of days in a year )
sea level | max number of freezing fronts [integer number]
Freezing front min and max depth [m]
saturation coefficient (fraction of 1)

grid.txt includes number of grid point (n), 
vertical grid (in the example starts from the 1.5 meter above the ground up to 90 meters) 
minus sign correspond to the values above the ground surface and pluus corresponds to the values below the surface.
Let me denote the n elements of the grid as "gr"
The n+1 element of the grid correspond to the number of output points (in the example it is 12)
The rest of the grid file correspond to indexes of the grid points (e.g. the number 40 below the 12 is the index of the 
gr(40)=0.001)

initial.txt includes all initial tempratures 
the code reads initial.txt file when in cmd.txt the first element is equal to 1. (e.i. start from the initial time step)
first parameter in the file can be ignored / second parameter is the number of points (in the example it is 13)
First column corresponds to the depth in meters and the second column to the temprature [Celsius] measured at that depth at time=0

geo.txt is the most important file which include all thermo-physical properties of the multilayered soil column
first row can be ignored
second row first element can be also ignored / second element in this row correspond to the number of layer 
starting from row 3 to row 8 are thermo-physical properties of each layer.
first column is the volumetric water content (WVC) [fraction of 1]
second and third columns are "a" and "b" coefficients of the unforzen water curve (obtained from unfrozen water curve fitting) [dimensionless]
forth and fifth columns are the thawed and frozen volumetric heat capacities [J/(m^3K)]
six and seven columns are the thawed and frozen heat conductivities [W/(mK)]
eigths column is the thickness of the corrsponding layer

Files input.txt and vegetaion.txt do not do anything. They needed for the code to run and was created for the parallel version of the code.


