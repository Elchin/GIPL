# GIPL (Geophysical Institute Permafrost Laboratory) Model

GIPL is a permafrost numerical model that employs phase changes and the effect of unfrozen volumetric water content in non-homogeneous soil texture. 
The original version of the model was originally developed by Romanovsky and Tipenko (2004) at the University of Alaska Fairbanks, and described in Marchenko et al., (2008). This version has been significantly modified from its predecessor and was adapted to the IRF (Initialize, Run, Finalize) coding standard structure (see [Basic Model Interface](http://csdms.colorado.edu/wiki/BMI_Description)). 
This version is maintained by [Elchin Jafarov](https://www.woodwellclimate.org/staff/elchin-jafarov/). Please cite Jafarov et al., (2012) when using the model. This model is a part of [permamodel project](https://github.com/permamodel/permamodel).

# References:
Jafarov, E. E., Marchenko, S. S., and Romanovsky, V. E.: Numerical modeling of permafrost dynamics in Alaska using a high spatial resolution dataset, The Cryosphere, 6, 613–624, https://doi.org/10.5194/tc-6-613-2012, 2012.

Jafarov, E. E., Romanovsky, V. E., Genet, H., McGuire A., D., Marchenko, S. S.: The effects of fire on the thermal stability of permafrost in lowland and upland black spruce forests of interior Alaska in a changing climate, Environmental Research Letters, 8, 035030, 2013. https://doi.org/10.1088/1748-9326/8/3/035030

Jafarov, E.E., Nicolsky, D.J., Romanovsky, V.E., Walsh, J.E., Panda, S.K., Serreze, M.C. 2014. The effect of snow: How to better model ground surface temperatures. Cold Regions Science and Technology, Volume 102, Pages 63-77, ISSN 0165-232X, [doi: 10.1016/j.coldregions.2014.02.007](http://www.sciencedirect.com/science/article/pii/S0165232X1400038X). 

Zlotnik, V.A., Harp, D.R., Jafarov, E.E., Abolt, C.J. A Model of Ice Wedge Polygon Drainage in Changing Arctic Terrain. Water 2020, 12, 3376. [https://doi.org/10.3390/w12123376](https://www.mdpi.com/2073-4441/12/12/3376)

# Installation: 

**Windows**: Compile the `gipl.f90` using [gfortran](https://gcc.gnu.org/wiki/GFortran) or similar compiler and name the executable file `gipl.exe`.

**Linux**: Use Makefile to create an executable. Navigate to the GIPL folde and type `make` in the terminal.

**Mac**:
```bash 
conda env create --file gipl.yml
source activate gipl
make
```
# Run: 
Make sure to create a folder called `out` before running the executable file (see [`gipl_config.cfg`](https://github.com/Elchin/GIPL/blob/master/gipl_config.cfg)).  <br />

# Visualize:
The file with measured temperatures is [`mesres.txt`](https://github.com/Elchin/GIPL/blob/master/mesres.txt). The header for the `mesres.txt` can be found in  the `compare.m` script. Type the command `>>compare(0)` to execute the script in the [Matlab](https://www.mathworks.com/products/matlab.html). The script plots the daily measured against calculated ground temperatures at four specified depth as shown in the figure below. Also checkout the jupyter notebook [example](https://github.com/Elchin/GIPL/blob/master/plot_results.ipynb).
![results](https://github.com/Elchin/GIPL/blob/master/results.png)

# Input/Output Setup:
[`gipl_config.cfg`](https://github.com/Elchin/GIPL/blob/master/gipl_config.cfg) configuration file that includes paths to input and output files. All input files stored in `in` folder. The path for the output files is prescribed in the config file. The model outputs three files: `results.txt` file with daily output, `mean.txt` with the yearly averaged data such as active layer thickness and freeze-up day, and the file `start.txt` includes the temperature profile with depth for the last day of simulation. The `result.txt` has the structure: `time`, `upper_bnd_temperature`, `snow_depth, ground_temperatures`, where 'ground_temperatures' assigned in [`grid.txt`](https://github.com/Elchin/GIPL/blob/master/in/grid.txt) (see below). The `mean.txt` file has the same configuration and includes 3 more columns. Everything in `mean.txt` is averaged yearly (see `number_of_time_steps` in `gipl_config.cfg`).

## Input data:
All input files are located in the [`in`](https://github.com/Elchin/GIPL/tree/master/in) folder.

### **gipl_config.cfg** 
Includes paths for input and output files and the correspoding setup paramteres. <br />
```
0/1: start from previous time step / start from the begining
step | taum | tmin :
    step is the timestep in the example it is 1
    taum is the convergence parameter used by the stefan subroutine 
    tmin minimal timestep used in the Stefan subroutine
begin | end : start and end, in the example it runs over one year from 0 to 1
smoothing_factor | unfrozen_water_parameter | max number of iterations
number_of_second_per_day [sec] | number_of_time_steps (in the example number of days in a year )
sea_level | max_number_of_freezing_fronts [integer number]
freezing_front_min | freezing_front_max depth [meters]
saturation_coefficient (fraction of 1)
```

### **input.txt** 
Includes the total number of sites and the corresponding ids for the organic and mineral soils. In the current version, the number of sites is equal to 1.

### **bound.txt** 
Include upper boundary condition (in the example it is an air temperature)<br />
```
First row is a number of observations (in the example number of day)
Column 1: Time parameter (day number)
Column 2: Temperature (daily averaged temperature [Celsius])
```

### **snow.txt** 
Include snow depth (in the example it is daily snow depth)<br />
```
First row is a number of observations (in the example number of day)
Column 1: Time parameter (day number)
Column 2: Snow depth (daily averaged [m])
```

### **rsnow.txt** 
Include snow thermal conductivity (in the example it is daily snow conductivity)<br />
```
First row is a number of observations (in the example number of day) 
Column 1: Time parameter (day number)
Column 2: Snow conductivity (daily averaged [W/(mK)])
```

### **grid.txt** 
Includes number of grid point (`n`), <br />
In the example, the vertical grid starts from the 1.5 meters above the ground up to 90 meters deep. 
The minus sign corresponds to the values above the ground surface and plus corresponds to the values below the surface.
For more clarity, copy and paste the grid into the excel file. The `n+1` element of the grid corresponds to the number of output points (in the example it is 12). The rest of the grid file correspond to indexes of the grid points (e.g. the number 40 below the 12 is the index of the gr(40)=0.001).

### **initial.txt** 
Includes the initial temperatures with depth profile.
The model reads `initial.txt` file when in `cmd.txt` the first element is equal to 1 (e.i. start from the initial time step).
The first parameter in `initial.txt` can be ignored, the second parameter is the number of points (in the example it is equal 13).
The first column corresponds to the depth [m] and the second column to the temperature [Celsius] measured at that depth at time=0.

### **mineral.txt** 
Includes all thermo-physical properties of the multilayered soil column.

The first row can be ignored.
In the second row, the first element can be also ignored. The second element in this row corresponds to the number of layers, 
starting from row 3 to row 8 are thermo-physical properties of each layer.
The first column is the volumetric water content (WVC)  is a fraction of 1.
The second and the third columns are "a" and "b" coefficients of the unfrozen water curve (obtained from unfrozen water curve fitting) [dimensionless].
The forth and the fifth columns are the thawed and frozen volumetric heat capacities [J/(m^3K)].
The six and the seven columns are thawed and frozen heat conductivities [W/(mK)].
The eighths column is the thickness of the corresponding layer.


### **organic.txt** 
Includes a similar structure to `mineral.txt` and carries the parameters for the organic soil layer/s.

### **sites.txt** 
This file includes `site_id`, `snow_code`, `veg_code`, `geo_code`, `gt_zone_code`, and `temp_grd`. This file is useful when we run multiple sites at the same time and have multiple organic and mineral soils, snow parameters, and so on. The current setup uses a temperature gradient equal 0.0 at the lower boundary of the grid.
