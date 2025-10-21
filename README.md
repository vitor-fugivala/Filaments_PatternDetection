Cosmic Web Reconstruction algorithm (CWR) developed by [Chen et al. (2015)](https://academic.oup.com/mnras/article/454/1/1140/1138949) for filament detection in a 2D galaxy distribution using the Subspace Constrained Mean Shift (SCMS) method.

The procedure is divided into two parts. First, we run **pattern_detection.py** to detect the filaments. Then, **bootstrap.py** is used to estimate their uncertainties through bootstrap resampling.

**_IMPORTANT: Bootstrap code runs CWR in parallel. Not to use more than 4 cores if RAM < 1 TB._**



# 1. Compiling core C++ scripts

First, we begin by compiling the main C++ code:

````
g++ kdtree_2d.hpp kdtree_2d.cpp find_filament_2d.cpp -o find_filament_2d
````
_Note: requires the [Eigen](https://github.com/PX4/eigen/tree/master) library._



# 2. Detecting cosmic filaments

Next, the **pattern_detection.py** script is executed with four arguments in the following order:

```
python3 pattern_detection.py <data> <bandwidth> <iter> <epsilon>
```

- `<data>:` Two-column `.txt` file containing the [RA, DEC] positions of all galaxies.
- `<bandwidth>:` _A_<sub>0</sub> parameter (`float`) defining the smoothing scale ([Eq. A1](https://academic.oup.com/mnras/article/454/1/1140/1138949)). Smaller values produce more filaments. Recommended range: [0.4, 0.8].
- `<iter>:` Maximum number of iterations (`int`).
- `<epsilon>:` ε threshold (`float`) — points stop moving once their displacement is less than this value. The algorithm terminates when all points have converged.

## Example input

```
python3 pattern_detection.py '/path/to/CWR/data.txt' 0.8 1000 1E-4
```

## Output files

- `realfil.txt:` [RA, DEC] coordinates of the detected filaments.
- `filaments.png:` Visualization of the filaments.



# 3. Bootstrap uncertainty estimation

Finally, the **bootstrap.py** script uses the output **realfil.txt** file to estimate uncertainties in the detected filaments. Be aware that this step is computationally expensive and may take several hours.

The same arguments described previously are used:

```
python3 bootstrap.py <data> <bandwidth> <iter> <epsilon>
```

There are two additional optional parameters:

- `-nj, --njobs:` number of CPU cores to use (`int`). Default: (total system cores) − 2
- `-nb, --nboots:` number of bootstrap samples (`int`). Default: 100

## Example input

```
python3 bootstrap.py '/path/to/CWR/data.txt' 0.8 1000 1E-4 -nj 1 -nb 50
```

## Output files

- `bootstrap_sample_i.txt:` resampled datasets, where _i_ ranges from 1 to _nboots_.
- `bootstrap_fil_i.txt:` detected filaments for the _i_-th bootstrap sample.
- `finalfil.csv:` file containing the filament coordinates and their associated errors.  
  Columns are, respectively: **RA**, **DEC**, **Error**, **Reject** (if *Reject = 1*, the detection is unstable and the point should be disregarded).  
- `finalfil.png:` visualization of the stable filaments with error bars.



# Note: running the standalone C++ code

If you wish to test the _find_filament_2d_ code independently, six parameters must be specified.

### Example input

```
./find_filament_2d -d '/path/to/CWR/data.txt' -h 0.230705 -i 10000 -t 0.0182435 -e 0.00001 -n filaments01.txt
```

- `-d:` Two-column `.txt` file containing the [RA, DEC] positions of all galaxies.
- `-h:` Smoothing bandwidth parameter (`float`).
- `-i:` Maximum number of iterations (`int`).
- `-t:` Threshold parameter (`float`) controlling the density level used for filament detection.
- `-e:` ε threshold (`float`).
- `-n:` output filename (`string`) for the detected filaments.
