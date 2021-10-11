SOFA - Spatially Oriented Format for Acoustics
==============================================

SOFA is a file format for reading, saving, and describing spatially
oriented data of acoustic systems.

Examples of data we consider are head-related transfer functions (HRTFs),
binaural room impulse responses (BRIRs), multichannel measurements such as done
with microphone arrays, or directionality data of loudspeakers.

The format specifications are the major issue, but we also aim in providing APIs
for reading and writing the data in SOFA.

For more information on the format specifications and available data have a look
at http://www.sofaconventions.org/


Downloads
=========

Current versions of SOFA can be found on its [old
home](http://sourceforge.net/projects/sofacoustics/files/?source=navbar).

At the moment we are working on a new release which will be shortly available on
this site.


Usage
=====

## Matlab/Octave API

In order to use SOFA with Matlab or Octave you have to add its `API_MO` folder
to your paths. After that you can play around with your acoustic measurements
as shown by the following example which uses a head-related transfer function
measurement.

```matlab
%% put your information here:
hrtf = SOFAload('path/to_your/HRTF.sofa');
soundInput = audioread('path/to_your/fancy_audio_file.wav');

%% demo script
% Start SOFA
SOFAstart;
% Display some information about the impulse response
SOFAinfo(hrtf);
% Plot a figure with the measurement setup
SOFAplotGeometry(hrtf);
% Have a look at the size of the data
disp(['size [MxRxN]: ' num2str(size(hrtf.Data.IR))])
% Calculate the source position from a listener point of view
apparentSourceVector = SOFAcalculateAPV(hrtf);
% Listen to the HRTF with azimuth of -90°
apparentSourceVector(91, 1)
SOFAplotGeometry(hrtf, 91);
soundOutput = [conv(squeeze(hrtf.Data.IR(91, 1, :)), soundInput) ...
               conv(squeeze(hrtf.Data.IR(91, 2, :)), soundInput)];
sound(soundOutput, hrtf.Data.SamplingRate);
```
