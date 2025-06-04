# Analysis FLow

This README explains the analysis logic flow. All variables are defined in the helper.py file.

## Step 1: miniAOD -> nano
Todo
## Step 2: nano -> flat ntuples
Todo
## Step 3: Skimming of flat ntuples
The used flat ntuples are indicated in the helper.py file as `<channel>_cons_25` and located at ```T3_PSI_CH```, under ```/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano```.
We skim these files using the skimmer located at: ```/RDsTools/skim```, by running:

```
conda activate tf #or any other python3 env with matching packages
create_skimmer.py <date_time> <channel> <selection> <cons> <prod>
```

where ```<date_time>``` corresponds to the date and time of the variable ```<channel>_cons_25``` in helper.py, ```<channel>``` is the corresponding channel, the ```<selection>``` must be defined in the ```baselines``` dictionary in the ```helper.py```, ```<cons>``` specifies constrained or unconstrained prouction (True for the final analysis), ```<prod>``` the production year (25 for the final analysis). Example:

```
create_skimmer.py 24_04_2025_14_03_22 data base_wout_tv_25 True 25
```
The skimmed files are located at ```/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/```, with the same datetime string as the unskimmed ones.

## Step 4: Calculate BDT weights for data

BDT weights are calcualted on data only to correct for mismodelings in the combinatorial background estimate. The BDT trainer is located at:
```/RDsTools/classification``` and is trained on the skimmed data files. It is run using:

```
#f.e.: book some worker nodes and run this interactively on the worker nodes
conda activate tf #need ML environment with XGBoost
python bdt_trainer.py -p <prod>
```

also here ```<prod>``` specifies the production year (25 for the final analysis). The used bdt model is specified in helper under ```bdt_data_25```. The BDT qualitiy and performance plots are stored at ```/RDsTools/classification/<bdt_data_25>```.
The bdt model can be evaluated using the evaluator stored in the same directory. It is run using:

```
python create_submitter_bdt.py -d <date_time> -s <sidebands> -p <prod> -n <nfiles> 

```

where ```<sidebands>``` specifies if we use the bdt trained on 'left' or both ('double') sidebands ('double' for final analysis). The BDT weighted data will be stored at ```/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data`/<date_time>```, with <date_time> being the ```<bdt_data_25>``` or ```/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data`/leftSB/<date_time>``` if evaluating only on the leftSB trained BDT.


## Step 5: Calculate Hammer weights for signals

For the hammer procedure we need unfiltered gen-level samples for all signals, and an additional dsmu simulation with the ISGW2 model for the circualar test. Similary to the other samples, they are indicated in the `helper.py` file as ```<channel>\_gen``` and saved under `/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano`.

Before calculating the hammer weights, we perform the circular test on the unfiltered gen-level dsmu signal samples. First, we calculate the weights using `hammer\_circular.cc` in `RDsTools/hammercpp/development\_branch\weights`:

```
conda activate hammercpp
cmake .
make
./submit_circular.sh <channel> <input> <target> <nfiles> <prod>
```
So run:
```
./submit_circular.sh dsmu       CLN   ISGW2 10000 25
./submit_circular.sh dsmu_isgw2 ISGW2 CLN   10000 25
```
The weighted samples are stored under ```/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/<prod>/``` and indicated in the `helper.py` file under: ```dsmu_isgw2_to_cln``` and ```dsmu_to_isgw2```. 
To produce the plots use `circularTest.cc` in `RDsTools/hammercpp/tests`:
```
./run_circularTest
```
The plots are saved under a datetime folder in the same directory.
Secondly, we now have to calculate the average unfiltered gen-level weights when weighting the unfiltered gen-level MC samples to BCL/BGL by using the `hammer\_temp.cc` in the same folder:

```
./submit_hammer.sh <channel> <input> <target> <nfiles> <prod>
```
So run:
```
./submit_hammer.sh dsmu       CLN   BCLVar 10000 25
./submit_hammer.sh dstau      ISGW2 BCLVar 10000 25
./submit_hammer.sh dsstarmu   CLN   BGLVar 10000 25
./submit_hammer.sh dsstartau  ISGW2 BGLVar 10000 25
```
The weighted samples are stored under ```/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/<prod>/``` and indicated in the `helper.py` file under: ```<channel>_to_<target>```. 
To calculate the average weights use `getAverageWeight.cc`:

```
`./run_averageWeight
```
The average weights will be stored as a yaml file in a datetime folde rin the same directory.

Finally, we weight our flat signal MC samples using again `hammer\_temp.cc`:


```
./submit_hammer.sh <channel> <input> <target> <nfiles> <prod>
```
So run with the default option, which reweights the different signals correctly to BCL/BGL (only posible for <channel> = "signal"):

```
./submit_hammer.sh signal default default 100000 25
```

The weighted samples are stored under ```/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/<prod>/``` and indicated in the `helper.py` file under: ```sig_cons_hammer_25```.
To produce the weight effect plots for the signal use `plotVariations.cc` in the same folder:

```
./run_plotVariations <var>

```
where ```<var>``` is the variable you want to plot. The plots are saved under a ```plots/<datetime>``` folder in the same directory.


