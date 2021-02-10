# compressive-sensing-Gerchberg-Saxton
Software implementation of multiple algorithms for the generation of 3D multi-foci patterns with phase-only spatial light modulators, for applications including optogenetics, neural imaging and optical tweezers.

While this library is good for a basic implementation and for understanding the algorithm, if you have access to a recent Nvidia Graphics card for your SLM setup, we strongly suggest using our GPU-based alternative to this library for much faster computation times (https://github.com/ppozzi/SLM-3dPointCloud).

The software is compatible with both Python 2.x and 3.x, and its only external dependency is numpy. The code consists in 5 functions implementing different algorithms for the calculation of holograms for three dimensional optogenetics stimulation with phase only spatial light modulators.

Running the Python file shows an example calculating computer generated holograms for a random three dimensional distribution of focal points in a 100x100x10 micromenters volume with a 20x microscope objective. More details and practical suggestions are available in the code comments.

The implemented algorithms are:

- Random superposition (RS)
- Gerchberg-Saxton (GS)
- Weighted Gerchberg-Saxton (WGS)
- Compressive Sensing Gerchberg-Saxton (CSGS)
- Weighted Compressive Sensing Gerchberg-Saxton (WCSGS)

RS, GS and WGS are implemented as reported in the paper by Di Leonardo et Al. (https://doi.org/10.1364/OE.15.001913).
A complete description of CSGS and WCSGS is presented in our publication "Fast Calculation of Computer Generated Holograms for 3D Photostimulation through Compressive-Sensing Gerchberg–Saxton Algorithm" (https://doi.org/10.3390/mps2010002), and on our upcoming publication "Real time generation of three dimensional patterns for multiphoton stimulation" (https://doi.org/10.3389/fncel.2021.609505 , https://www.frontiersin.org/articles/10.3389/fncel.2021.609505/full).

If you use this code for academic research, please consider citing the reported literature.
