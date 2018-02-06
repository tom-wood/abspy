# Abspy
A series of macros for TOPAS and python scripts for adjusting for the 
effects of absorption in X-ray diffraction samples using a capillary in 
Debye-Scherrer geometry. The relevant literature is <i>J. Appl. Crystallogr.
</i>, 2015, <b>48</b>, 149; 2017, <b>50</b>, 240 (and references therein).

The macros are based around a "brute-force" method of high-order polynomial 
fitting (hence the large number of coefficients).

## Getting Started
Adjusting for the intensities is straightforward: all that is needed is the
macro\_abs\_intensities.txt file, which contains the relevant TOPAS macro.
Although this macro (and the other macros) are written as if you can refine
the value of muR, in practice this needs to be done manually (with TOPAS 6
at least).

Adjusting for 2 theta offsets is also very straightforward, but only if the
capillary used had a diameter of 1 mm or 0.7 mm. In these cases, the macros
already exist as macro\_th2offset\_1.txt or macro\_th2offset\_pt7.txt and
can be pasted into the TOPAS input file. Obviously the value of muR should
be the same as in the intensities macro. Also the secondary radius (Rs) of 
the diffractometer is required.

To get the macros, just click on the files and then right-click "Raw" and
"Save link as".

## More complicated things
### Creating a bespoke th2offset macro
This involves running a python script and so requires a working python
installation along with the following libraries: numpy
