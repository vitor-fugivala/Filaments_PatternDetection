# Filaments_PatternDetection

Cosmic Web Reconstruction algorithm (Chen et al. 2015) for filament detection in a 2D galaxy distribution using the Subspace Constrained Mean Shift (SCMS) method.

The procedure is divided into two parts. First, we run pattern_detection.py to detect the filaments. Then, bootstrap.py is used to estimate their uncertainties through bootstrap resampling.


1)
Scripts necessaries to execute the Cosmic Web Reconstrucion code (Chen et al. 2015).
The scripst are in C++. The main one is in Python ("pattern_detection.py"), that connects all the others. 
IMPORTANT: Pyhton code runs CWR in parallel. Not to use more then 4 cores if RAM < 1TB.

#-------------------------------------------------------------------

We begin by compiling the main C++ code:

Bash

g++ kdtree_2d.hpp kdtree_2d.cpp find_filament_2d.cpp -o find_filament_2d
Note: requires the Eigen library.

#-------------------------------------------------------------------

Next, the pattern_detection.py script is run with 4 arguments:
```
python3 pattern_detection.py -d -h -i -e
```
-d: .txt file with two columns containing the [RA, DEC] positions of all galaxies.

-h: A_0 parameter for calculating the smoothing/bandwidth parameter h (Eq. A1 [https://arxiv.org/abs/1501.05303]; the paper recommends values in the range [0.4, 0.8]).

-i: Maximum number of iterations.

-e: Epsilon — points stop moving once the movement is less than epsilon. The algorithm terminates if all points stop moving.

INPUT EXAMPLE for the pattern_detection.py script:
```
python3 pattern_detection.py '/path/to/CWR/data.txt' 0.8 1000 1E-4
```
The script returns:

realfil.txt — [RA, DEC] coordinates of the detected filaments

filaments.png — visualization of the filaments

#-------------------------------------------------------------------

Finally, the bootstrap.py script uses the realfil.txt file to calculate uncertainties. By default, the data is resampled 100 times. This step is computationally expensive and may take several hours.

The same arguments used previously must be provided.
```
python3 bootstrap.py -d -h -i -e
```

INPUT EXAMPLE for the bootstrap.py script:
```
python3 bootstrap.py '/path/to/CWR/data.txt' 0.8 1000 1E-4
```

The script returns:

bootstrap_sample_i.txt — resampled datasets, where i is in the range [1, 100].

bootstrap_fil_i.txt — detected filaments for the bootstrap_sample_i.txt sample.

finalfil.csv — final file containing the filaments and their errors. The columns are, respectively: RA, DEC, Error (uncertainty determined by the bootstrap method), err_threshold (maximum accepted error value), Reject (if Reject == 1, the point should be disregarded. If Reject == 0, the detection is stable).

finalfil.png — visualization of the filaments with errors.

#-------------------------------------------------------------------

Note: to run the standalone C++ program:

INPUT EXAMPLE for find_filament_2d:
```
./find_filament_2d -d '/path/to/CWR/data.txt' -h 0.230705 -i 10000 -t 0.0182435 -e 0.00001 -n filaments01.txt
```
(Using h and t values calculated by the "pattern_detection.py" script. The program also includes the -n argument for the output filename)
