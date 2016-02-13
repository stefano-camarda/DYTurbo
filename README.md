# DYTURBO - fast Drell-Yan Monte Carlo and quadrature integrator
DYTURBO is a fast Drell-Yan integrator, based on the DYRES and DYNNLO programs.
Numerical integration is performed using both the Vegas Monte Carlo integration and quadrature rules.

## How to setup USER
1. ask for guest rights and get tarbal from webpage
```
https://gitlab.cern.ch/DYdevel/DYTURBO/wikis/home
tar xzvf dyturbo-VERSION.tar.gz
cd dyturbo-VERSION
```

2. setup and compile
```
(optional PATH=<lhapdfdir>:$PATH if you want to use your version of lhapdf previously installed in <lhapdfdir>)
(optional setup root to enable root output)
./configure --enable-root --enable-Ofast
make
make install
```

3. have fun
```
./bin/dyturbo
```

## How to setup DEVEL
 1. ask for developer access and checkout the repository
```
git clone https://gitlab.cern.ch/DYdevel/DYTURBO.git
```

 2. setup/compile
```
mkdir m4
autoreconf -i
(optional PATH=<lhapdfdir>:$PATH if you want to use your version of lhapdf previously installed in <lhapdfdir>)
./configure
make && make install
```

 3. have fun (by default with CT10NLO)
```
./bin/dyturbo
```

## Description of project
 - autotools are used for building project
 - on first configuration `LHAPDF` and `Cuba` will be downloaded and installed
     (if not found already)

 - additional configuration options:
    - `--enable-O3`             Use compiler optimization flags -O3.
    - `--enable-Ofast`          Use unsafe compiler optimization flags -Ofast.
    - `--enable-root`           Use root for histograming. (default=no)

 - debug options:
    - `--enable-debug`          Add debug flags to include gdb debug symbols.
    - `--enable-trapFPE`        Stop of floating point errors.
    - `--enable-checkBounds`    add -fbounds-check flag for compilation.


## Directory structure
- `cernlib` contains code from CERNLIB libraries
- `dynnlo` contains code from DYNNLO
- `dyres`  contains code from DYRES
- `mcfm`  contains code from MCFM
- `src` contains DYTURBO C++ code which is steering the calculations and
     running the Fortran procedures
- `scripts` contains scripts for plotting and submitting to 
- input files can be found in `input` folder (where else :) )


## Customization
 - two files can be change to customize output of the calculation:
     - `src/plotter.C` for histograms
     - `src/settings.C` function `cuts` for changing cuts

## Merging results
 1. Merging of jobs
     - we recommend to run DYTURBO per each term (RES,CT,FO,REAL,VIRT) separately
     - to obtain correct normalization from run of several jobs with different
       random seeds, please, use our merger
    ```
    ./bin/merger -X merged_file.root result_file_1.root result_file_2.root ...
    ```
     - Note: in case of real term for NNLO prediction it is better to use median results
       (object with suffix median in output file of `merger` program). It will remove ouliers
       from your distribution, which are caused by color dipole cancellations
 2. Merging of terms
     - if you want to obtain finite order prediction you should sum up counter
       term you should sum up fixed order (i.e. FO for NLO and REAL+VIRT for NNLO)
       to counter term (simply with hadd)
     - to obtain final prediction you need to finite order to resummation part ( simply with hadd)
