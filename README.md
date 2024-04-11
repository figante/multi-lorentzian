Here you'll find an .xcm file for Xspec [example_constant_phase_lags,xcm](https://github.com/figante/multi-lorentzian/blob/main/example_constant_phase_lags.xcm) to fit simultaneously:    

Power spectrum band 1 (PS1)    
Power spectrum band 2 (PS2)    
Real part of the cross spectrum of band 1 vs. band 2 (RE)    
Imaginary part of the cross spectrum of band 1 vs. band 2 (IM)    

The .xcm file uses a linear combination of Lorentizans to fit PS1 and PS2    
and a combination of Lorentzians multiplied by the cosine (sine) of a parameter    
(the phase lag) to RE (IM).    

The .xcm file also reads the phase-lag (PL) spectrum and the coherence function (CF)    
and computes the derived model of these two on the basis of the model fitted to   
PS1, PS2, RE and IM.    

See for [Mendez et al. 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.9405M/abstract) for details    
