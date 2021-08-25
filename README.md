# steerablepyrtexture: Steerable pyramid decomposition and texture analysis tools 

A modification pyrtools which is the Python 3 port of Eero Simoncelli's matlabPyrTools.
This version only includes the steerable pyramid portions of pyrtools and adds the texture
analysis function from their textureSynth MATLAB package. Additionally, to make 


## Dependencies

From sptools:
 - numpy
 - scipy
 - matplotlib
 - Pillow
 - tqdm
 - requests



# Authors

Original Authors:

Rob Young and Eero Simoncelli, 7/13

William Broderick, 6/17

William Broderick, Pierre-Étienne Fiquet, Zhuo Wang, Zahra Kadkhodaie,
Nikhil Parthasarathy, and the Lab for Computational Vision, 4/19

Modified by:
Nicholas Rogers 6/21

# Usage:

To generate a texture feature vector from a 256x256 image the command

```
import steerablepyrtexture as spt
params = spt.textureAnalyze(img, 5, 4, 7, vector_output=True)
```

will use a 5 subband, 4 orientation steerable pyramid decomposition with a 7x7 square 
of entries around the center of the auto/cross-correlation matrices. 

The default output is a dictionary containing the parameter descriptors. Use vector output 
as a shortcut for generating a feature vector. With the settings shown above the output will 
have 2195 features. Some of which will should be always 0.




