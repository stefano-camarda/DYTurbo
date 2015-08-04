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
./configure
make && make install
```

 3. have fun (by default with CT10NLO)
```
./bin/dyfast
```

## Description
 - input files can be found in `input` folder (where else :) )


