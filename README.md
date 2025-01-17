# szm
Here is a code for astronomical source image extraction procedures (C++) and its application to meteor event detection using Python 3 (Astropy + NumPy + SciPy). 

Detection of meteor event images using the SSA400 telescope at the Assy-Turgen Observatory

I.Izmailov(1,2), M.Khovrichev(1,2), A.Tolstoy(2), S. Pavlov(2), D. Bikulova(1,2), M. Krugov(3), S. Sittykova(3)

(1) The Central Astronomical Observatory of the RAS at Pulkovo 
(2) Institute of Applied Astronomy of the Russian Academy of Sciences 
(3) Fesenkov Astrophysical Institute, Aerospace Committee of the Ministry of Digital Development, Innovations and Aerospace Industry of the Republic of Kazakhstan


The analysis of space mission data, such as LDEF, reveals that most meteor particles fall within the microgram to milligram range. This range corresponds to meteor event absolute magnitudes from 5-6 mag to 10-11 mag. The data processing based on the results of the most successful meteor monitoring networks (GMN, MMT, GWAC) demonstrates incompleteness of the detected event number in this magnitude range. Hence, developing a telescope with a 40-50cm aperture, a focal ratio of approximately 1:1, and a field of view of several tens of square degrees is a natural progression for the evolution of meteor event monitoring systems. Real-time data processing is necessary for telescopes and cameras due to the significant excess of data compared to typical storage volumes. Therefore, improving the meteor detection algorithms is relevant for separate images and a set of frames. This paper describes our method for detecting meteor events using the SSA400 system (D/F = 400/551 mm). The telescope was developed by astronomers of the Fesenkov Astrophysical Institute of the Aerospace Committee of the Ministry of Digital Development, Innovations and Aerospace Industry of the Republic of Kazakhstan and installed in the Assy-Turgen Observatory. For now, the telescope is used in a wide range of astronomical programs, including monitoring meteor events. The proposed algorithm consists of three stages: 1) Extracting all possible source images in each frame. 2) Applying the Hough transform to this set of sources to detect meteor tracks. 3) Performing peak detection on the set of diagrams in Hough's space to confirm the detection. The paper (https://doi.org/10.31725/0367-7966-2024-234-22-35, https://ui.adsabs.harvard.edu/abs/2024PPulO.234...22I/abstract) proves the high efficiency of this approach.

Installation and usage.

1. Copy all files including the fits directory into a folder, make the folder a current (cd folder), and compile a lib (libszm.so) by: 

g++ -fPIC -shared -o libszm.so szm.c

or

g++ -std=gnu++17 -fPIC -shared -o libszm.so szm.c

2. Download sample files via http://outreach.puldb.ru/fits.zip , and unzip these files to the fits folder.

Next, there are two examples of usage...

3. Run the ht-and-detect-in-frame-set.py by:

python3 ht-and-detect-in-frame-set.py

This will show a file name that contains a meteor event and create a plot with the detection peak.

4. Run ht-and-detect.py by:

python3 ht-and-detect.py

This will provide meteor line parameters and save several plots for visualization

Both examples are presented as a *.ipynb file for the jupyter-notebook.
