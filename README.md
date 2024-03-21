<h1><img src="logo/qhptomography.png" alt="logo" width='215' align="right"/> ebeDREENA</h1>

ebeDREENA is computational framework for generating high-pT predictions based on a dynamical energy loss formalism. The framework can include any, in principle arbitrary, event-by-event fluctuating temperature evolution within the dynamical energy loss formalism. This version is generalized to account for both LHC and RHIC energies and collision systems.

## < 1 > compilation

Compilation of the source code, performed using gcc compiler:

g++ source/*.cpp -fopenmp -O3 -o ebeDREENA

a) -fopenmp is necessary to enable parallelization using OpenMP;  
b) -O3 optimization is not necessary, but recommended;  
c) -o DREENAA is optional and if omitted, the output of the compilation will be placed in a.out;  

## < 2 > prerequisite files

All prerequisite files need to be textual tables or binary files. Textual tables can have a different number of columns depending on the file and they can contain header lines that start with '#' symbol;  

#### a) initial pT distributions

initial pT distributions file should have 2 columns in format:

pT | dsigma/d(pT^2) |
--- | --- |
... | ... |

for heavy flavor, initial pT distributions can be obtained from this [web site](http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html);  
> [!CAUTION]
> note that the default output of this web interface is dsigma/dpT, while ebeDREENA initial pT distribution input needs to be dsigma/d(pT^2), so these distributions need to be modified;

to avoid modifying different parameters that are hard-coded, initial pT distribution for heavy flavour should be in the range of 1GeV to at least 200GeV and for light flavour from 1GeV to at least 450GeV; for most distributions pT step of 1GeV seems to be sufficient;  

initial pT distribution file path relative to the executable should be: `./pTDists/ptDists[sNN]/ptDist_[sNN]_[particleName].dat`, where *sNN* is collision energy that can be 200GeV, 2760GeV, 5020GeV and *particleName* is the name of the particle that can be Bottom, Charm, Down, DownBar, Gluon, Strange, Up, UpBar;  

#### b) binary collision points

binary collision points path relative to the executable should be: `./binarycollpts/binarycollpts_cent=[centrality]/binarycollpts[eventID].dat`, where *centrality* is centrality parameter in format 'xx-xx%' (ie 0-5%, 10-20%,...), and *eventID* is non negative integer that represents event's ID;


binary collision points are generated by initial condition model and are used to generate jet's initial positions in x-y plane, while direction angles are uniformly distributed;


within this repository there is an example of binary collision points file generated using Monte Carlo Glauber initial condition model (see [arxiv:2208.09886](https://inspirehep.net/literature/2139859) for more details);  

#### c) temperature evolutions

temperature evolution files should be in binary format as 32bit floating point numbers;  

temperature evolution file path relative to the executable should be: `./evols/evols_cent=[centrality]/tempevol[eventID].dat`, where *centrality* is centrality parameter in format 'xx-xx%' (ie 0-5%, 10-20%,...), and *eventID* is non negative integer that represents event's ID;

within this repository there is an example temperature evolution file generated evolving Monte Carlo Glauber initial condition (also provided in this repository) with 3D hydro for 5020GeV Pb+Pb collision (see [arxiv:2208.09886](https://inspirehep.net/literature/2139859) for more details);

#### d) temperature evolution grid parameters file

temperature evolution grid parameters file is located in the same directory as the temperature evolution files: `./evols/evols_cent=[centrality]/evolgridparams.dat`, where *centrality* is centrality parameter in format 'xx-xx%' (ie 0-5%, 10-20%,...);

this file contains temperature evolution grid parameters sush as tau0, x and y grid ranges and steps; this format for evolutions is chosen to save storage space and import time; within this repository there is an example temperature evolution grid parameter file that corresponds to provided temperature evolution;

### e) LTables

> [!CAUTION]
> note that the default output of this web interface is dsigma/dpT, while ebeDREENA initial pT distribution input needs to be dsigma/d(pT^2), so these distributions need to be modified;  

to avoid modifying different parameters that are hard-coded, initial pT distribution for heavy flavour should be in the range of 1GeV to at least 200GeV and for light flavour from 1GeV to at least 450GeV; for most distributions pT step of 1GeV seems to be sufficient;  
initial pT distribution file path relative to the executable should be: ./pTDists/ptDists[sNN]/ptDist_[sNN]_[particleName].dat, where [sNN] is collision energy that can be 200GeV, 2760GeV, 5020GeV and [particleName] is the name of the particle that can be Bottom, Charm, Down, DownBar, Gluon, Strange, Up, UpBar; 

#### c) binary collision density

binary collision density file contains jet creation probability in transverse plane as a function of x and y; x and y need to form ordered grid with run order y > x, while probability is the 3rd column in the table;  
the binary collision density can be given only in the 1st quadrant of the transverse plane (x>=0, y>=0) if initial conditions have this symmetry;  
binary collision density needs to be consistent with initial conditions used to generate temperature evolution;  
binary collision density provided in this demo (./binarycolldensities/binarycolldensity_cent=30-40%.dat) is generated using optical Glauber model (see [arxiv:2110.01544](https://inspirehep.net/literature/2606181) for more details) and it corresponds to provided temperature evolution;

#### b) temperature evolution

temperature evolution file contains temperature as a function of proper time, tau, and x and y spatial coordinates in the transverse plane in that order; tau, x and y need to form an ordered grid with run order y > x > tau;  
temperature evolution file can contain an additional column with energy density evolution; if the file contains temperature and energy density evolution, the order of these two columns can be arbitrary (energy density can be 4th column and temperature 5th, or vice versa);  
the evolution can be given only in the 1st quadrant of the transverse plane (x>=0, y>=0) if initial conditions have this symmetry;  
the evolution provided in this demo (./evols/tempevol_cent=30-40%.dat) is generated evolvong optical Glauber initial conditions using 3D hydro model (see [arxiv:2110.01544](https://inspirehep.net/literature/2606181) for more details);

#### d) LTables

LTables files contain pre-generated radiated gluon rates and collisional energy loss; for radiative energy loss, there are 2 tables lnorm table with column format:

tau | p | temp | LNorm
--- | --- | --- | --- |
... | ... | ... | ... |

and ldndx table with column format:

tau | p | temp | x | Ldndx
--- | --- | --- | --- | --- |
... | ... | ... | ... | ... |

these tables are also a function of, effective number of flavours, (nf=2.5 for sNN=200GeV and nf=3.0 for LHC collision energies), particle mass (particle name is in file name) and chromo-magnetic to electric mass ratio, xB, which figure in file names;

for collisional energy loss, there is one table lcoll with column format:
p | temp | LColl
--- | --- | --- |
... | ... | ... |

this table is also a function of particle mass (particle name is in file name) and effective number of flavours, nf;

LTables files path relative to the executable should be:

+ `./ltables/lnorm_nf=[nf]_[particleName]_xB=[xB].dat`,
+ `./ltables/ldndx_nf=[nf]_[particleName]_xB=[xB].dat` and
+ `./ltables/lcoll_nf=[nf]_[particleName].dat`,

where *nf* is the effective number of flavours that can be 2.5 for 200GeV or 3.0 for LHC energy collisions, *particleName* is the name of the particle that can be Bottom, Charm, LQuarks or Gluon (all light quarks are taken to have same mass, so their LTables are the same), and *xB* is the chromo-magnetic to electric mass ratio;

unlike previous files, ebeDREENA calculates LTables; however, these tables need to be calculated only once and can be reused while calculating high-pT energy loss with different temperature evolution backgrounds;

within this repository there is an example of LTables files for Charm quark, nf=3.0 and xB=0.6;  

#### d) phiGausPts

in `./phiGaussPts/` directory are textual tables containing jet's direction angles and weights that correspond to Gaussian quadrature integration method in range [0, 2Pi];

jet's direction angles are sampled in these points, so that afterwards, when integrating final pT,phi distribution over phi, which is nedeed to obtain R_AA, v_2, v_3,..., there is no need for angle resampling;

## < 3 > running ebeDREENA  

There are two possible calculation options within ebeDREENA framework: LTables calculation and energy loss calculation. Since ebeDREENA is parallelized using OpenMP, set OMP_NUM_THREADS environmental variable to desired value before running calculations.

#### a) LTables calculation

to see all parameters and their default values run:

```
./ebeDREENA LTables -h
```

parameters for LTables calculations are:

+ **sNN** parameter: case sensitive string with possible options: 200GeV, 2760GeV, 5020GeV, 5440GeV;  
*default value: PbPb*

+ **pName** parameter: case sensitive string with possible options: Bottom, Charm, LQuarks, Gluon, where LQuarks stand for light quarks, since all light quarks (down, down-bar, strange, strange-bar, up and up-bar) use the same LTables;  
*default value: Charm*

+ **xB** parameter: float that represents magnetic to electric mass ratio; based on lattice calculation: xB=0.6;  
*default value: 0.6*

additional parameters are number of points used for QuasiMonteCarlo integration LdndxMaxPoints and LCollMaxPoints, with default values 500000 and 10000, respectively;  

to generate LTables provided in this demo, that are for 5020GeV collision energy, charm quark and for xB value of 0.6, use:

```
./ebeDREENA LTables --sNN=5020GeV --pName=Charm --xB=0.6

```

or just:

```
./ebeDREENA LTables

```

since these are all default parameter values;

#### b) energy loss calculation

to see all parameters and their default values run:  

```
./ebeDREENA AverageEL -h
```

parameters for energy loss calculations are:  

+ **collsys** parameter: case sensitive string with possible options: AuAu, PbPb, XeXe,...  
*default value: PbPb*

+ **sNN** parameter: case sensitive string with possible options: 200GeV, 2760GeV, 5020GeV, 5440GeV;  
*default value: 5020GeV*

+ **pName** parameter: case sensitive string with possible options: Bottom, Charm, Gluon, Down, DownBar, Strange, Up, UpBar, LQuarks;  
the calculation can be done for each parton individualy, however there is an option to calculate all light quarks at the same time using modified algorithm since the only thing differentiating light quarks is initial pT distribution; this leads to 4x speed-up of the calculation time compared to calculating each light quark individually;  
*default value: Charm*

+ **centrality** parameter: string in format 'xx-xx%' (ie 0-5%, 10-20%,...);  
*default value: 30-40%*

+ **xB** parameter: float representing magnetic to electric mass ratio; based on latest lattice calculation: xB=0.6;  
*default value: 0.6*

+ **eventN** parameter: non-negative integer representing number of events;  
number of binary collision points and temperature evolutions should match this number; event loop goes from 0 up to not including *eventN*:  
`for (size_t eventID=0; eventID<eventN; ++eventID)`;  
*default value: 1000*

+ **BCPP** parameter: string representing binary collision points percentage, for example '25%';  
this parameter determines the percentage of binary collision points that will be used as jet's initial positions, meaning that events with larger number of binary collisions will have more jet's since the percentage is the same;  
*default value: 20%*

+ **phiGridN** parameter: non-negative integer representing number of jet's direction angles;  
for each initial position, *phiGridN* jets are generated;  
based on *phiGridN* parameter, angle values are imported from `./phiGaussPts/phiptsgauss[phiGridN].dat`  
*default value: 25*

+ **TIMESTEP** parameter: positive float, that represents the timestep of jet traversing qgp medium, expressed in fm;  
*default value: 0.1*

+ **TCRIT** parameter: positive float, that represents critical temperature expressed in GeV, i.e. the temperature value for which the energy loss stops;  
*default value: 0.155*

+ **BCPSEED** parameter: non-negative integer that represents seed for random engine that samples binary collision points as jet's initial positions;  
when set to 0, no seed is set and every run will produce different results;  
*default value: 0*

## < 4 > outputs of ebeDREENA

Output of LTables run is in ./ltables/ directory; for single run, there are three output files - two for radiative and one for collisional. These tables are determined by effective number of flavours, i.e. collision energy, particle type and the value of xB. These tables are necessary for energy loss calculations and they are reused for different different events, centralitites hydro backgrounds obtained with different models,... As an example, tables for charm quark are given for LHC collision energy and for xB value of 0.6. 

Output of energy loss run is in `./results/results[particleName]/`, where *particleName* can be Bottom, Charm, Down, DownBar, Gluon, Strange, Up, UpBar. These are R_AA(pT,phi) distributions that are later used to calculate R_AA and v_n. These files also have headers, that contain various informations about the event such as event ID, average jet path-length and temperature jets experience along the trajectory,...

For file name pattern see file provided in this repository.  

R_AA(pT,phi) distribution provided in `./results/resultsCharm/` is generated using initial pT distribution obtained from http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html, as well as binary collision points, temperature evolution and ltables provided in this repository with all default parameter values except *BCPSEED*, which is set to 1, so that users can check if they get the same results.  

## < 5 > disclaimer

This version of the ebeDREENA framework is meant for high-pT energy loss calculation on event-by-event fluctuating hydro background.  
ebeDREENA framework provides results for bare quarks and gluons. To obtain results comparable to experimental data, fragmentation functions such as DSS, BCFY and KLP have to be used.  
Calculation time for LTables calculation is large (up to about hour and a half on 112 cores), however, they need to be calculated only once. Energy loss calculatio time for single event on single core is up to about 3 minutes, however calculation time vary for different particles.  

---------------------------------------------------------------
for aaditional questions contact Zigic Dusan at zigic@ipb.ac.rs
