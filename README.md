SOFA - Spatially Oriented Format for Acoustics
==============================================

SOFA is a file format for reading, saving, and describing spatially
oriented data of acoustic systems.

Examples of data we consider are head-related transfer functions (HRTFs),
binaural room impulse responses (BRIRs), multichannel measurements such as done
with microphone arrays, or directionality data of loudspeakers.

The format specification is the major focus of SOFA, but we also aim in providing 
toolboxes for reading and writing the data in SOFA. For more information on the 
format specifications and available data, see http://www.sofaconventions.org/.

This project implements a reference toolbox for SOFA: The SOFA Toolbox. 

SOFA Toolbox 2.x implements SOFA versions 2.x. The SOFA Toolbox has been previously 
known as the SOFA API_MO, which supported SOFA versions up to 1.x.


Downloads
=========

Current releases of SOFA Toolbox can be found on its [old
home](http://sourceforge.net/projects/sofacoustics/files/?source=navbar).


Usage
=====

## Matlab/Octave 

In order to use SOFA with Matlab or Octave add its `SOFAtoolbox` folder
to your search paths. After that you can play around with your acoustic measurements
as shown by the following example which uses a HRTF measurement.

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
