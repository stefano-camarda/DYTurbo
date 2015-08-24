# DYTURBO - fast Drell-Yan Monte-Carlo
DYTURBO is fast Drell-Yan Monte-Carlo based on DYRES code. Several improvements
has been done to increase the speed of calculation.

## How to setup
 1. checkout the repository
```
git clone https://gitlab.cern.ch/DYdevel/DYTURBO.git
```

 2. setup/compile
```
autoreconf -i
(optional PATH=<lhapdfdir>:$PATH if you want to use your version of lhapdf previously installed in <lhapdfdir>)
./configure
make && make install
```

 3. have fun (by default with CT10NLO)
```
./bin/dyfast
```

## Description of project
 - autotools are used for building project
 - on first configuration `LHAPDF` and `Cuba` will be downloaded and installed
     (if not found already)
 - additional option to configuration:
     - `--enable-debug` to compile with GDB debug symbols
     - `--enable-root` to save outputs in ROOT histograms
     - `--enable-fast`  use compiler optimization flags -O3

- `dyres` contains optimised Fortran code
 - `src` contains DYTURBO C++ code which is steering the calculations and
     running the Fortran procedures
 - `scripts` contains scripts for plotting and submitting to 
 - input files can be found in `input` folder (where else :) )


## Customization
 - two files can be change to customize output of the calculation:
     - `src/plotter.C` for histograms
     - `src/settings.C` function `cuts` for changing cuts


