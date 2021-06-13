# ContinuousPhaseModulation (CPM) - EuclideanDistance

A Matlab code for the calculation of the minimum squared Euclidean distance (d_min^2), based on the algorithm found in the Book: __Digital Phase Modulation__ page 463-464, for different types of CPM signals (GMSK-RECT-RC....).


# How to cite this Work ?
If you use this package for your own work, please consider citing it with [this piece of BibTeX](PERvsSPEC.bib):

```bibtex
@ARTICLE{9406053,
  author={Kassan, Karim and Farès, Haïfa and Christian Glattli, D. and Louët, Yves},
  journal={IEEE Transactions on Communications}, 
  title={Performance vs. Spectral Properties For Single-Sideband Continuous Phase Modulation}, 
  year={2021},
  volume={},
  number={},
  pages={1-1},
  doi={10.1109/TCOMM.2021.3073792}
  }
```



# How to run the code ?
1. Make sure that you have a compatibel version of matlab (this code was tested on matlab 2018b)
2. Download (clone) the files from the repository.
3. Open the file called _CPM_EuclideanDistance.m_
4. Select the section called __Pulse shape & Variable ini__
5. Select the type of pulses by changing the variable `pulse` number
	* `1` is for Lorentzian
	* `2` is for GMSK
	* `3` is for Raised Cosine
	* `4` is for Rectangular
6. Change the pulse length by changing the variable `L`
7. Select the sampling frequency by changing the variable `Fs` (Usually 64 is high enough)
8. Select M-ary by changing the `M` ( e.g M=2 for Binary)
9. Increase the variable `h_max` to increase the modulation index domain (h) in the plot figure.
10. Select the number of observation symbols (N)  ( increasing N-> dmin appraoch the upper bound)
	* Go to section __Main code__ to __Minimum Euclidien Distance__ part.
	* Select the variable name `Nmax`



# Example

In this part we reproduce the figure 3.23 page 93 of the Book : __Digital Phase Modulation__ for observation symbols $N=5$
````
pulse         = 4;
L             = 1;
Fs            = 64;
M             = 2^1;
h_max         = 1.2;
Nmax          = 5;
````
Results Plot:
![](Figures/EuclideanDistance_Rec_N5.png?raw=true)
## License
© 2018-2019 Karim Kasan

