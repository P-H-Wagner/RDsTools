# Create Tensorflow GPU env to run with NVIDIA-SMI 495.29.05    Driver Version: 495.29.05    CUDA Version: 11.5 

```
conda create -n gpu python=3.9 -y
conda activate gpu

conda install -c conda-forge cudatoolkit=11.2 cudnn=8.2 -y
pip install tensorflow==2.10
conda install pandas matplotlib scikit-learn seaborn pyyaml -c conda-forge
conda install -c conda-forge root=6.24.06
conda install numpy=1.23
pip install root_numpy

```
Before using the gpu environment, run this (or better, put in .sh file in activate.d folder):
```
export PYTHONNOUSERSITE=1 #to pick numpy from gpu env
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH #to pick cuda from gpu env
```

# To create a XGBoost GPU env instead, exchange ```pip install tensorflow==2.10``` with: 
```
pip install xgboost==1.6.1 --extra-index-url https://pypi.nvidia.com
```
 
