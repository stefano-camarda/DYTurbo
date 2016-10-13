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
git clone ssh://git@gitlab.cern.ch:7999/DYdevel/DYTURBO.git
git co devel
```

 2. setup/compile and test
```
autoreconf -i
(optional PATH=<lhapdfdir>:$PATH if you want to use your version of lhapdf previously installed in <lhapdfdir>)
./configure --enable-test --enable-Ofast [--enable-root]
make && make install && make check
```

 3. branch your repository add your code and dont forget to test ;)

        make check


## Description of project
 - autotools are used for building project
 - on first configuration `LHAPDF` and `Cuba` will be downloaded and installed
     (if not found already)

 - additional configuration options:
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
     - to obtain final prediction you need to finite order to resummation part (simply with hadd)


## Merge program description
 - To print basic usage just run with help argument

         [user@host DYTURBO]$ ./bin/merger -h
         usage: -h [-h] [-v] [-x]  <output> <input list>
          Please keep separated switches!!! Its on my todolist! 
             -x    Normalize histograms to Xsection. 
             -v    Increase verbosity (very chatty). 
             -h    Print this help message. 

          For more info read README.md
          ./bin/merger --help
          usage: --help [-h] [-v] [-x]  <output> <input list>
           Please keep separated switches!!! Its on my todolist! 
              -x    Normalize histograms to Xsection. 
              -v    Increase verbosity (very chatty). 
              -2d   Make 2d projections and outliers for 2D. 
              -h    Print this help message. 

          For more info read README.md


### Detailed description

- basic naming convention: object could be `TH1x`, `TH2x`, `TProfile`, `TProfile2D`. When
  I say average or median I mean per each bin.

- merger program produces several stuff on one merge:
    - for `TH1x`: average (no suffix), median (suffix `_median`)

    - for `TH1x` with name `h_qt` and `pt`: it creates also rebinned version
      (suffix `_rebin`) with Zpt measurement binning and their mean and median
      version

    - for `TH2x`: mean and median same as TH1x + creates projection to each axis
      (`_px`, `_py`) and their mean and median version

    - for `TProfiles`: average, i.e. `TProfile::Add(o,1/N)` of all object (no
      suffix) and outlier removal (suffix "outlier"). More details bellow.

    - for `TProfiles2D`: average + my plan is to add projections so we can fill only 2D during dyturbo run.

- There is option `-X`, which turns on "Scaling to final cross section": At the
  end of dyturbo I am saving result of integration into extra 2D histogram per
  each term and qt-y bin. This information is then used to correct integral of
  `TH1x`, `TH2x`.

### Outlier removal for profiles
- Right now, outliers are removed for `TProfile`s and preparing also for 2D
  (for other objects there is median version). 
- each bin of profile consist of two terms: nominator `sum(w*y)` ( the first moment `y` of
  cros section) and denominator `sum(w)`  (cross section), where `y` is studied observable,
  `w` is weight of event and we are suming per all events.
- in first step we create denominator median `sumhat(w_i) = median ( sum(w_ij)`
  from all merged files `j`, per each bin `i`
- we calculate chi-square from denominator

                    (sum(w_ij)-sumhat(w_i))**2
          chi2_ij = --------------------------
                     (delta (sum(w_ij)) )**2

- then we decide about outlier by `Prob(chi2_ij,1) > 7 sigma ` throw away
- from the rest we get files with indexes `l`
- final `outlier` profile values terms are: 

          sum(w_i*y_i) = sum(w_il*y_il) / N_i
          sum(w_i) = sum(w_il) / N_i

  where `N_i` is number of jobs contributing to bin `i`.
- the `TProfile` is then created with hard-set of above values


## Submission script
Basic script for submission is `submit_DYTURBO.sh`. It could be customize to
submit jobs to any batch system. If you need help, please run:

     ./scripts/submit_DYTURBO.sh --help

Or contact me [Jakub.Cuth@cern.ch]. An example for grid is here

    ./scripts/submit_DYTURBO.sh --grid --proc wp,wm,z0 --collider lhc7 --order 2 --term RES,CT,REAL,VIRT --pdfvar all

Dont forget to configure and compile before submiting

    ./configure --enable-Ofast --enable-root
    make -j
    make install


There are several grid related scripts in `scripts` folder, which you can edit for your need 

* `run_prun.sh`: function for submit jobs with prun command
* `compile_grid.sh`: compilation on grid side
* `run_grid.sh`: running on grid side


