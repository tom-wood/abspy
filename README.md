# Abspy
A series of macros for TOPAS and python scripts for adjusting for the 
effects of absorption in X-ray diffraction samples using a capillary in 
Debye-Scherrer geometry. The relevant literature is <i>J. Appl. Crystallogr.
</i>, 2015, <b>48</b>, 149; 2017, <b>50</b>, 240 (and references therein).

The macros are based around a "brute-force" method of high-order polynomial 
fitting (hence the large number of coefficients).

## Getting Started
Adjusting for the intensities is straightforward: all that is needed is the
[macro\_abs\_intensities.txt](../master/macro_abs_intensities.txt) 
file, which contains the relevant TOPAS macro.
Although this macro (and the other macros) are written as if you can refine
the value of muR, in practice this needs to be done manually (with TOPAS 6
at least).

Adjusting for 2 theta offsets is also very straightforward, but only if the
capillary used had a diameter of 1 mm or 0.7 mm. In these cases, the macros
already exist as [macro\_th2offset\_1.txt](../master/macro_th2offset_1.txt) or [macro\_th2offset\_pt7.txt](../master/macro_th2offset_pt7.txt) 
and can be pasted into the TOPAS input file. Obviously the value of muR 
should be the same as in the intensities macro. Also the secondary radius 
(Rs) of the diffractometer is required.

To get the macros, just click on the files and then right-click "Raw" and
"Save link as".

## More complicated things
### Creating a bespoke th2offset macro
This involves running a python script and so requires a working python
installation along with the following libraries: numpy, scipy, pandas, 
numba and matplotlib (anaconda or canopy are possible distributions).

The two python scripts required are [abs\_peak\_pos.py](../master/abs_peak_pos.py) and [abs\_th\_offset.py](../master/abs_th_offset.py). 
With the former, the values that should be changed are
the r value (radius in mm) and the save\_fname value. Running the script
will then output a series of profiles (saved as csv files) for each 2
theta value. The abs\_th\_offset.py script can then be used (once rad,
macro\_fname and fnames have been set) to output a macro with the
pertinent values for the capillary radius in question.

### Other files
The other files are added for completeness and are for those who want to 
check/modify the methods used. They comprise:
1. [abs\_fit.py](../master/abs_fit.py), which attempts to provide a 
user defined convolution (thus doing both intensity and 2 theta offset in
one), however, didn't work very well (very slow and was buggy in TOPAS).
2. [muR\_intensities.txt](../master/muR_intensities.txt), which is 
taken from the first reference.
3. [muR\_intensities.py](../master/muR_intensities.py), which fits
those intensities to a series of polynomials, which form the basis of the
intensity correction macro.
4. [absorption\_intensities.txt](../master/absorption_intensities.txt) 
is a version of the intensity correction macro, but not in macro form (was
tested using an IFDEF statement).
5. [absorption\_intensities\_rewrite.py](../master/absorption_intensities_rewrite.py) is a python script to rewrite 4 into a macro.
