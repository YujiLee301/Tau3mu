# Tau3mu
This repository privides the approach to generate Tau3mu pure Phase Space process with Pythia8.
## Setup ##
1. download Pythia8 and hepMC3 package.
2. install hepMC3 firstly, Reference: https://gitlab.cern.ch/hepmc/HepMC3
```
wget http://hepmc.web.cern.ch/hepmc/releases/HepMC3-3.2.6.tar.gz
  tar -xzf HepMC3-3.2.6.tar.gz
  mkdir hepmc3-build
  cd hepmc3-build
  cmake -DCMAKE_INSTALL_PREFIX=../hepmc3-install   \
        -DHEPMC3_ENABLE_ROOTIO:BOOL=OFF            \
        -DHEPMC3_ENABLE_PROTOBUFIO:BOOL=OFF        \
        -DHEPMC3_ENABLE_TEST:BOOL=OFF              \
        -DHEPMC3_INSTALL_INTERFACES:BOOL=ON        \
        -DHEPMC3_BUILD_STATIC_LIBS:BOOL=OFF        \
        -DHEPMC3_BUILD_DOCS:BOOL=OFF     \
        -DHEPMC3_ENABLE_PYTHON:BOOL=ON   \
        -DHEPMC3_PYTHON_VERSIONS=2.7     \
        -DHEPMC3_Python_SITEARCH27=../hepmc3-install/lib/python2.7/site-packages \
        ../HepMC3-3.2.6
  make
  make install
```
3. Compile Pythia with hepMC as the external Package
```
cd pythia8310
./configure --with-hepmc3=/mnt/e/hepmc3-install
make
```
4. Compile Tau3mu.cc with g++
```
g++ Tau3mu.cc -o Tau3mu -I../include -O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread  -L../lib -Wl,-rpath,../lib -lpythia8 -ldl -I/mnt/e/hepmc3-install/include -L/mnt/e/hepmc3-install/lib -Wl,-rpath,/mnt/e/hepmc3-install/lib -lHepMC3 -DHEPMC3
```
5. Interactive Running
```
./Tau3mu
```
6. Get fastsim root file with Delphes(Please install Delphes firstly)
```
./DelphesHepMC3 cards/delphes_card_CMS.tcl ../pythia8310/Delphes_CMS_Tau3mu.root ../pythia8310/examples/hepmcouttau.dat
```
