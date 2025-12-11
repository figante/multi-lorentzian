Here you'll find six spectra (with corresponding matrices) and two .xcm files for Xspec to fit simultaneously:

Power spectrum band 1 (PS1)    
Power spectrum band 2 (PS2)    
Real part of the cross spectrum of band 1 vs. band 2 (RE)    
Imaginary part of the cross spectrum of band 1 vs. band 2 (IM)    

The .xcm file for the [constant phase-lags model](https://github.com/figante/multi-lorentzian/blob/main/example_constant_phase_lags.xcm) uses a linear combination of Lorentzians 
to fit PS1 and PS2 and a combination of Lorentzians multiplied by the cosine (sine) of a 
parameter (the phase lag) to fit the RE (IM) part of the cross spectrum.    

The .xcm file for the [constant time-lags model](https://github.com/figante/multi-lorentzian/blob/main/example_constant_time_lags.xcm) uses a linear combination of Lorentzians 
to fit PS1 and PS2 and a combination of Lorentzians multiplied by the cosine (sine) of a 
parameter (the time lag) times 2 \pi the frequency (in Xspec the energy) to fit the RE (IM) part of the cross spectrum.    

Both .xcm files also read the phase-lag (PL) spectrum and the coherence function (CF)    
and compute the derived model of these two on the basis of the model fitted to   
PS1, PS2, RE and IM.    

The .xcm files provided have many comments explaining what they do at each step

See for [Mendez et al. 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.9405M/abstract) for details    

Download the .tar.gz file with all the data and .xcm files for the example
