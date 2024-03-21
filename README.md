# multi-lorentzian
# .xcm file for Xspec to fit simultaneously:
#
# Power spectrum band 1 (PS1)
# Power spectrum band 2 (PS2)
# Real part of the cross spectrum of band 1 vs. band 2 (RE)
# Imaginary part of the cross spectrum of band 1 vs. band 2 (IM)
# The .xcm file uses a linear combination of Lorentizans to fit PS1 and PS2
# and a combination of Lorentzians multiplied by the cosine (sine) of a parameter
# (the phase lag) to RE (IM).
# The .xcm file also reads the phase-lag (PL) spectrum and the coherence function (CF)
# and computes the derived model of these two on the basis of the model fitted to 
# PS1, PS2, RE and IM.
#
# See for [Mendez et al. 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.9405M/abstract) for details

statistic chi

# Here I read the data. I read 6 datasets: Numbers 1 and 2 are PS1 and PS2 
# (in the soft and hard bands, respectively), 3 and 4 are, respectively
# the RE and IM of the cross spectrum using the same energy bands as for PS1
# and PS2, and 5 and 6 are, respectively, PL and CF for the same 2 bands.

# One thing to notice is the way I define the responses. By doing it
# that way I can define, below, 6 different models, one for each of the
# 6 different datasets. Read the Xspec manual for this (named models,
# data command and response command)

# You will have to use the names of your own files here

data 1:1 pds_band1.pha
response  1:1 pds_band1.rmf

data 1:2 pds_band2.pha
response  1:2 none
response  2:2 pds_band2.rmf

data 1:3 real_band1-band2.pha
response  1:3 none
response  3:3 real_band1-band2.rmf

data 1:4 imag_band1-band2.pha
response  1:4 none
response  4:4 imag_band1-band2.rmf

data 1:5 lags_band1-band2.pha
response  1:5 none
response  5:5 lags_band1-band2.rmf

data 1:6 cohe_band1-band2.pha
response  1:6 none
response  6:6 cohe_band1-band2.rmf

# use here your own ignore command
# ignore 1:1-3,293-432 2:1-3,293-432 3:1-3,293-432 4:1-3,293-432 5:1-3,293-432 6:1-3,293-432

# These are standard Xspec settings.
method leven 10 0.01
abund angr
xsect vern
cosmo 70 0 0.73
xset delta 0.01
xset LINECRITLEVEL  1e-12

# This is the Lorentzian for the Real part, as in eq. 7 of my paper

mdefine xlor Lorentz(LineE,Width)cos(plag) : add

# This is the Lorentzian for the Imaginary part, as in eq. 7 of my paper

mdefine ylor Lorentz(LineE,Width)sin(plag) : add

# This is the definition of the coherence as in eq. 10 of my paper

# Because the way it is defined, I need to put inside more
# Lorentzians than the ones I will use. This is so because if I want to add another one
# later, I cannot use "addcomponent" in xspec, and I have to define
# the whole model again.
# By doing it this way I can put the norms of all the unecessary Lorentzians to 0 at the start, and freeze all the parameters of the Lorentzian, and start letting them free if I need to
# add a new Lorentzian.

mdefine coherence ((xlor(LineE1,Width1,plag1)*norm1+xlor(LineE2,Width2,plag2)*norm2+xlor(LineE3,Width3,plag3)*norm3+xlor(LineE4,Width4,plag4)*norm4+xlor(LineE5,Width5,plag5)*norm5+xlor(LineE6,Width6,plag6)*norm6+xlor(LineE7,Width7,plag7)*norm7+xlor(LineE8,Width8,plag8)*norm8+xlor(LineE9,Width9,plag9)*norm9+xlor(LineE10,Width10,plag10)*norm10)**2+(ylor(LineE1,Width1,plag1)*norm1+ylor(LineE2,Width2,plag2)*norm2+ylor(LineE3,Width3,plag3)*norm3+ylor(LineE4,Width4,plag4)*norm4+ylor(LineE5,Width5,plag5)*norm5+ylor(LineE6,Width6,plag6)*norm6+ylor(LineE7,Width7,plag7)*norm7+ylor(LineE8,Width8,plag8)*norm8+ylor(LineE9,Width9,plag9)*norm9+ylor(LineE10,Width10,plag10)*norm10)**2)/((Lorentz(LineE1,Width1)*nn1+Lorentz(LineE2,Width2)*nn2+Lorentz(LineE3,Width3)*nn3+Lorentz(LineE4,Width4)*nn4+Lorentz(LineE5,Width5)*nn5+Lorentz(LineE6,Width6)*nn6+Lorentz(LineE7,Width7)*nn7+Lorentz(LineE8,Width8)*nn8+Lorentz(LineE9,Width9)*nn9+Lorentz(LineE10,Width10)*nn10)*(Lorentz(LineE1,Width1)*mm1+Lorentz(LineE2,Width2)*mm2+Lorentz(LineE3,Width3)*mm3+Lorentz(LineE4,Width4)*mm4+Lorentz(LineE5,Width5)*mm5+Lorentz(LineE6,Width6)*mm6+Lorentz(LineE7,Width7)*mm7+Lorentz(LineE8,Width8)*mm8+Lorentz(LineE9,Width9)*mm9+Lorentz(LineE10,Width10)*mm10)) : add

# This is the definition of the lags as in eq. 9 of my paper
# The same considerations about the number of Lorentzians applies here as it applied in the coherence.

mdefine plags atan2(ylor(LineE1,Width1,plag1)*norm1+ylor(LineE2,Width2,plag2)*norm2+ylor(LineE3,Width3,plag3)*norm3+ylor(LineE4,Width4,plag4)*norm4+ylor(LineE5,Width5,plag5)*norm5+ylor(LineE6,Width6,plag6)*norm6+ylor(LineE7,Width7,plag7)*norm7+ylor(LineE8,Width8,plag8)*norm8+ylor(LineE9,Width9,plag9)*norm9+ylor(LineE10,Width10,plag10)*norm10,xlor(LineE1,Width1,plag1)*norm1+xlor(LineE2,Width2,plag2)*norm2+xlor(LineE3,Width3,plag3)*norm3+xlor(LineE4,Width4,plag4)*norm4+xlor(LineE5,Width5,plag5)*norm5+xlor(LineE6,Width6,plag6)*norm6+xlor(LineE7,Width7,plag7)*norm7+xlor(LineE8,Width8,plag8)*norm8+xlor(LineE9,Width9,plag9)*norm9+xlor(LineE10,Width10,plag10)*norm10) : add

# Here I define the model to fit the PDS in band1. This would be like eq. 3 (first row) in my paper.
# Notice that I only have 3 Lorentzxians. All the others are set to 0 and frozen.
# Also, I do not use 0, but 1e-10 as 0, because this makes the fits and the errors more stable. It
# doees not make any difference because 1e-10 is effectively 0, but the fits are more stable. 

model  1:1pds lorentz + lorentz + lorentz + lorentz + lorentz + lorentz + lorentz + lorentz + lorentz + lorentz
      0.0580861       0.01      1e-10      1e-10      1e+10      1e+10
      0.0377136       0.01      1e-10      1e-10      1e+10      1e+10
     0.00201022       0.01      1e-10      1e-10      1e+10      1e+10
        5.84269       0.01      1e-10      1e-10      1e+10      1e+10
        12.1826       0.01      1e-10      1e-10      1e+10      1e+10
     0.00903052       0.01      1e-10      1e-10      1e+10      1e+10
         4.4505       0.01      1e-10      1e-10      1e+10      1e+10
        7.57986       0.01      1e-10      1e-10      1e+10      1e+10
     0.00270649       0.01      1e-10      1e-10      1e+10      1e+10
      0.0167851         -1      1e-10      1e-10      1e+10      1e+10
      0.0470719         -1      1e-10      1e-10      1e+10      1e+10
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
       0.573699         -1      1e-10      1e-10      1e+10      1e+10
        1.84671         -1      1e-10      1e-10      1e+10      1e+10
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
    1.16505e-10         -1      1e-10      1e-10      1e+10      1e+10
    1.06076e-10         -1      1e-10      1e-10      1e+10      1e+10
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
        4.81115         -1      1e-10      1e-10      1e+10      1e+10
        3.19062         -1      1e-10      1e-10      1e+10      1e+10
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
       0.590984         -1      1e-10      1e-10      1e+10      1e+10
        1.76853         -1      1e-10      1e-10      1e+10      1e+10
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
      0.0106788         -1      1e-10      1e-10      1e+10      1e+10
      0.0253539         -1      1e-10      1e-10      1e+10      1e+10
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
        4.76518         -1      1e-10      1e-10      1e+10      1e+10
        3.79816         -1      1e-10      1e-10      1e+10      1e+10
          1e-10         -1      1e-10      1e-10      1e+10      1e+10

# Here I define the model to fit the PDS in band1. This would be like eq. 3 (second row) in my paper.
# Notice how some of the parameters are linked to PDS1
# I only leave the normalisations free

model  2:2pds lorentz + lorentz + lorentz + lorentz + lorentz + lorentz + lorentz + lorentz + lorentz + lorentz
= 1pds:p1
= 1pds:p2
      0.0100210       0.01      1e-10      1e-10      1e+10      1e+10
= 1pds:p4
= 1pds:p5
      0.0114238       0.01      1e-10      1e-10      1e+10      1e+10
= 1pds:p7
= 1pds:p8
     0.00058553       0.01      1e-10      1e-10      1e+10      1e+10
= 1pds:p10
= 1pds:p11
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
= 1pds:p13
= 1pds:p14
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
= 1pds:p16
= 1pds:p17
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
= 1pds:p19
= 1pds:p20
          1e-10         -1      1e-10      1e-10      1e+10      1e+100
= 1pds:p22
= 1pds:p23
          1e-10         -1      1e-10      1e-10      1e+10      1e+100
= 1pds:p25
= 1pds:p26
          1e-10         -1      1e-10      1e-10      1e+10      1e+100
= 1pds:p28
= 1pds:p29
          1e-10         -1      1e-10      1e-10      1e+10      1e+10

# Here I define the model to fit the Real part of the cross spectrum (first row of eq. 8)
# Notice how some of the parameters are linked to PDS1

model  3:3rea xlor + xlor + xlor + xlor + xlor + xlor + xlor + xlor + xlor + xlor
= 1pds:p1
= 1pds:p2
        1.75855       0.01      -1000      -1000       1000       1000
    0.000182521       0.01      1e-10      1e-10      1e+10      1e+10
= 1pds:p4
= 1pds:p5
       0.169791       0.01      -1000      -1000       1000       1000
      0.0115519       0.01      1e-10      1e-10      1e+10      1e+10
= 1pds:p7
= 1pds:p8
       -2.85097       0.01      -1000      -1000       1000       1000
     0.00348894       0.01      1e-10      1e-10      1e+10      1e+10
= 1pds:p10
= 1pds:p11
              0         -1      -1000      -1000       1000       1000
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
= 1pds:p13
= 1pds:p14
              0         -1      -1000      -1000       1000       1000
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
= 1pds:p16
= 1pds:p17
       -0.04245         -1      -1000      -1000       1000       1000
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
= 1pds:p19
= 1pds:p20
              0         -1      -1000      -1000       1000       1000
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
= 1pds:p22
= 1pds:p23
              0         -1      -1000      -1000       1000       1000
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
= 1pds:p25
= 1pds:p26
              0         -1      -1000      -1000       1000       1000
          1e-10         -1      1e-10      1e-10      1e+10      1e+10
= 1pds:p28
= 1pds:p29
              0         -1      -1000      -1000       1000       1000
          1e-10         -1      1e-10      1e-10      1e+10      1e+10

# Here I define the model to fit the Imaginary part of the cross spectrum (second row of eq. 8).
# Notice that there is no free parameters cause these are the same parameters of the Real part (see eq. 8), 

model  4:4ima ylor + ylor + ylor + ylor + ylor + ylor + ylor + ylor + ylor + ylor
= 3rea:p1
= 3rea:p2
= 3rea:p3
= 3rea:p4
= 3rea:p5
= 3rea:p6
= 3rea:p7
= 3rea:p8
= 3rea:p9
= 3rea:p10
= 3rea:p11
= 3rea:p12
= 3rea:p13
= 3rea:p14
= 3rea:p15
= 3rea:p16
= 3rea:p17
= 3rea:p18
= 3rea:p19
= 3rea:p20
= 3rea:p21
= 3rea:p22
= 3rea:p23
= 3rea:p24
= 3rea:p25
= 3rea:p26
= 3rea:p27
= 3rea:p28
= 3rea:p29
= 3rea:p30
= 3rea:p31
= 3rea:p32
= 3rea:p33
= 3rea:p34
= 3rea:p35
= 3rea:p36
= 3rea:p37
= 3rea:p38
= 3rea:p39
= 3rea:p40

# Here I define the model of the lags. Notice that there is no free parameter, because the values of the parameters come fromt the fit of PDS1 PDS RE and IM

model  5:5pla plags
= 3rea:p1
= 3rea:p2
= 3rea:p3
= 3rea:p4
= 3rea:p5
= 3rea:p6
= 3rea:p7
= 3rea:p8
= 3rea:p9
= 3rea:p10
= 3rea:p11
= 3rea:p12
= 3rea:p13
= 3rea:p14
= 3rea:p15
= 3rea:p16
= 3rea:p17
= 3rea:p18
= 3rea:p19
= 3rea:p20
= 3rea:p21
= 3rea:p22
= 3rea:p23
= 3rea:p24
= 3rea:p25
= 3rea:p26
= 3rea:p27
= 3rea:p28
= 3rea:p29
= 3rea:p30
= 3rea:p31
= 3rea:p32
= 3rea:p33
= 3rea:p34
= 3rea:p35
= 3rea:p36
= 3rea:p37
= 3rea:p38
= 3rea:p39
= 3rea:p40
              1         -1          0          0      1e+20      1e+24

# Here I define the model of the coherence. Notice that there is no free parameter, because the values of the parameters come fromt the fit of PDS1 PDS RE and IM

model  6:6coh coherence
= 3rea:p1
= 3rea:p2
= 3rea:p3
= 3rea:p4
= 3rea:p5
= 3rea:p6
= 3rea:p7
= 3rea:p8
= 3rea:p9
= 3rea:p10
= 3rea:p11
= 3rea:p12
= 3rea:p13
= 3rea:p14
= 3rea:p15
= 3rea:p16
= 3rea:p17
= 3rea:p18
= 3rea:p19
= 3rea:p20
= 3rea:p21
= 3rea:p22
= 3rea:p23
= 3rea:p24
= 3rea:p25
= 3rea:p26
= 3rea:p27
= 3rea:p28
= 3rea:p29
= 3rea:p30
= 3rea:p31
= 3rea:p32
= 3rea:p33
= 3rea:p34
= 3rea:p35
= 3rea:p36
= 3rea:p37
= 3rea:p38
= 3rea:p39
= 3rea:p40
= 1pds:p3
= 1pds:p6
= 1pds:p9
= 1pds:p12
= 1pds:p15
= 1pds:p18
= 1pds:p21
= 1pds:p24
= 1pds:p27
= 2pds:p30
= 2pds:p3
= 2pds:p6
= 2pds:p9
= 2pds:p12
= 2pds:p15
= 2pds:p18
= 2pds:p21
= 2pds:p24
= 2pds:p27
= 2pds:p30
              1         -1          0          0      1e+20      1e+24
bayes off

# Here I ignore dataeerts 5 and 6 because I am not fitting them. Obce I have the fit I can notice them again and the model of the lags and the coherence should pass through the data.

ignore 5-6:**

