Cosmic Web Reconstruction algorithm (CWR) developed by [Chen et al. (2015)](https://academic.oup.com/mnras/article/454/1/1140/1138949) for filament detection in a 2D galaxy distribution using the Subspace Constrained Mean Shift (SCMS) method.

The procedure is divided into two parts. First, we run **pattern_detection.py** to detect the filaments. Then, **bootstrap.py** is used to estimate their uncertainties through bootstrap resampling.

**_IMPORTANT: Bootstrap code runs CWR in parallel. Not to use more then 4 cores if RAM < 1TB._**



# 1. Compiling core C++ scripts

First, we begin by compiling the main C++ code:

````
g++ kdtree_2d.hpp kdtree_2d.cpp find_filament_2d.cpp -o find_filament_2d
````
_Note: requires the [Eigen](https://github.com/PX4/eigen/tree/master) library._



# 2. Detecting cosmic filaments

Next, the **pattern_detection.py** script is executed with four arguments in the following order:

```
python3 pattern_detection.py -d -h -i -e
```

- `-d:`  two-column `.txt` file containing the [RA, DEC] positions of all galaxies.
- `-h:` _A_<sub>0</sub> parameter ([Eq. A1](https://academic.oup.com/mnras/article/454/1/1140/1138949)) that defines the smoothing scale. Smaller values produce more filaments. Recommended range: [0.4, 0.8].
- `-i:` maximum number of iterations.
- `-e:` ε threshold. Points stop moving once their displacement is less than this value. The algorithm terminates when all points have converged.

## Example input:

```
python3 pattern_detection.py '/path/to/CWR/data.txt' 0.8 1000 1E-4
```

## Output files:

- `realfil.txt:` [RA, DEC] coordinates of the detected filaments
- `filaments.png:` visualization of the filaments



# 3. Bootstrap uncertainty estimation

Finally, the **bootstrap.py** script uses the **realfil.txt** file to calculate uncertainties. This step is computationally expensive and may take several hours.

The same arguments described previously must be provided:

```
python3 bootstrap.py -d -h -i -e
```

There are two optional parameters:

- `-nj, --njobs:` number of cores. Default: (total system cores) - 2
- `-nb, --nboots:` number of bootstrap samples. Default: 100

## Example input:

```
python3 bootstrap.py '/path/to/CWR/data.txt' 0.8 1000 1E-4
```

## Output files:

- `bootstrap_sample_i.txt:` resampled datasets, where i is in the range [1, 100].
- `bootstrap_fil_i.txt:` detected filaments for the bootstrap_sample_i.txt sample.
- `finalfil.csv:` — final file containing the filaments and their errors. The columns are, respectively: RA, DEC, Error (uncertainty determined by the bootstrap method), err_threshold (maximum accepted error value), Reject (if Reject == 1, the point should be disregarded. If Reject == 0, the detection is stable).
- `finalfil.png:` — visualization of the filaments with errors.



# Note: standalone C++ code testing

## Example imput

```
./find_filament_2d -d '/path/to/CWR/data.txt' -h 0.230705 -i 10000 -t 0.0182435 -e 0.00001 -n filaments01.txt
```

- `n:` output filename
