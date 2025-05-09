# Hammer setup

To use the hammer library (https://gitlab.com/mpapucci/Hammer) with C++, follow these steps.

1. Install conda 
2. Create a conda environment with the following packages (versions are important to avoid clashes later!)

```
conda create -n hammerpp root=6.22.6
conda activate hammercpp
conda install cmake=3.20.2
conda install boost=1.76.0 #Dont pick a newer boost version than this!
conda install -c conda-forge eigen
```

3. Set up hammer itself. We keep all examples, force the yamlcpp and hepmc installation as we don't have them, and install hammer with root. Choose static libraries ```(-DBUILD_SHARED_LIBS=OFF)``` over shared ones:

```
wget https://hammer.physics.lbl.gov/Hammer-1.4.1-Source.tar.gz
tar -xzf Hammer-1.4.1-Source.tar.gz

mkdir Hammer-build
cd Hammer-build

cmake -DCMAKE_INSTALL_PREFIX=../Hammer-install -DWITH_EXAMPLES=ON -DWITH_EXAMPLES_EXTRA=ON -DWITH_ROOT=ON -DFORCE_HEPMC_INSTALL=ON -DFORCE_YAMLCPP_INSTALL=ON -DENABLE_TESTS=ON -DINSTALL_EXTERNAL_DEPENDENCIES=ON -DBUILD_SHARED_LIBS=OFF ../Hammer-1.4.1-Source
```
Now we can run make, make install and perform an optinal test if wanted:

```
make 
ctest -V      # test the installation (optional)
make install
```

## Change from CLN to HQET2 (used in EvtGen)

The correspondance between Hammers CLN and EvtGens HQET2 parameters is defined in ```getCLNParameters``` in ```hammer_all.cc```. Since EvtGen sets ```f0m``` to zero (corresponds to ```Fz``` in Hammer), we have to modify the hammer source code in https://gitlab.com/mpapucci/Hammer/-/blob/master/src/FormFactors/CLN/FFBtoDCLN.cc?ref_type=heads#L104:

```
const double Fz = 0.0; 
```


