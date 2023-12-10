# js8call-lsb

Alberto Garlassi I4NZX December 2023
This code allows to receive the JS8 subband in LSB mode instead of USB.
Receiver should be tuned + 3 kHz from usual, eg. 7081 kHz for 40 m.
It's real use is for DSB receivers like ADX. Unfortunately an ADX tuned at 7078
for normal USB decoding also feeds JS8CALL with the mirrored FT8 subband, which starts at 7074.
Tuning 3 kHz up may bring in SSB traffic, but the power spectrum there is usually way less hot then for FT8
The RECEIVE_LSB directive is located in src/lib/js8/js8_params.f90

This code works only for reception. Transmission is standard USB and RIT should be set at +3kHz.
In other words: if you just want to check RX tune +3, for transmitting check if your RIT
has a 3 kHz span, tune the rig to the official frequency and set RIT at 3 kHz.
If RIT span is not enough, change the source code of your rig. Or maybe the split frequency option
of JS8Call could be used.
- 
The code for JS8Call was cloned on December 10, 2023 from the official BitBucket repo
https://bitbucket.org/widefido/js8call/src/js8call/

Only two files were changed:
src/lib/js8/js8_params.f90
src/lib/decoder.f90

Download the original 2.2.1-devel code, substitute the two files, follow the build instructions.
Tested on XUbuntu only.

If someone finds this code useful, please commit it to the main repo, maybe one day there will be a new JS8Call release.

Credits and copyright: please refer to the JS8Call page.