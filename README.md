# compressive-sensing-Gerchberg-Saxton
Software implementation of multiple algorithms for the generation of 3D optogenetics photostimulation patterns with phase-only spatial light modulators

The software is written for Python 2.7.x , and its only external dependency is numpy. The code consists in 5 functions implementing different algorithms for the calculation of holograms for three dimensional optogenetics stimulation with phase only spatial light modulators. Running the Python file shows an example calculating computer generated holograms for a random three dimensional distribution of focal points in a 100x100x10 micromenters volume with a 20x microscope objective. More details and practical suggestions are available in the code comments.

The implemented algorithms are:

- Random superposition (RS)
- Gerchberg-Saxton (GS)
- Weighted Gerchberg-Saxton (WGS)
- Compressive Sensing Gerchberg-Saxton (CSGS)
- Weighted Compressive Sensing Gerchberg-Saxton (WCSGS)

RS, GS and WGS are implemented as reported in the paper:
"Di Leonardo et Al. "Computer generation of optimal holograms for optical trap arrays." Optics Express 15.4 (2007): 1913-1922."
A complete description of CSGS and WCSGS will be presented in an upcoming publication.
