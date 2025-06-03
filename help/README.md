# Analysis FLow

This README explains the analysis logic flow. All variables are defined in the helper.py file.

## Step 1: miniAOD -> nano
Todo
## Step 2: nano -> flat ntuples
Todo
## Step 3: Skimming of flat ntuples
The used flat ntuples are indicated in the helper.py file as ```<channel>_cons_25````and located at ```T3_PSI_CH```, under ```/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano```.
We skim these files using the skimmer located at: ```/RDsTools/skim```, by running:

```
conda activate tf #or any other python3 env with matching packages
create_skimmer.py <date_time> <channel> <selection> <cons> <prod>
```

where ```<date_time>``` corresponds to the date and time of the variable ```<channel>_cons_25``` in helper.py, ```<channel>``` is the corresponding channel, the ```<selection>``` must be defined in the ````baselines``` dictionary in the ```helper.py```, ```<cons>``` specifies constrained or unconstrained prouction (True for the final analysis), ```<prod>``` the production year (25 for the final analysis). Example:

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

where ```<sidebands>``` specifies if we use the bdt trained on 'left' or both ('double') sidebands ('double' for final analysis). The BDT weighted data will be stored at ```/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data`/<date_time>```, with <date_time> being the ```<bdt_data_25>```.
