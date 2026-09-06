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

This project took advantage of netCDF software developed by NSF Unidata (http://doi.org/10.5065/D6H70CW6).

Downloads
=========

Current releases of SOFA Toolbox can be found on [SONICOM Ecosystem](https://ecosystem.sonicom.eu/tools/14).

**Note:** 
- The SOFA Toolbox 2.1 is the succeeding release version of SOFA API M/O version 1.1.3. It supports SOFA 2.1 as known as AES69-2022.
- The SOFA Toolbox 2.0 has never been released. 


Usage
=====

## Matlab/Octave 

In order to use SOFA with Matlab or Octave add its `SOFAtoolbox` folder
to your search paths. After that you can play around with your acoustic measurements
as shown by the following example which uses a HRTF measurement.

```matlab
%% Start SOFA Toolbox
% Add SOFA Toolbox directory to Matlab/Octave paths before starting
SOFAstart;

%% Load SOFA file
% The following command loads SOFA file 'dtf b_nh5.sofa' from local directory - if available - otherwise downloads it from the online SOFA database
Obj = SOFAload('db://database/ari/hrtf b_nh5.sofa');

%% Inspect the data
disp('HRTFs loaded - Inspect the data structure and continue with "dbcont" when done.');
keyboard;

%% Display some information about the SOFA file
% Show general meta data, and measurement info
SOFAinfo(Obj);

%% Plot geometry of measurement setup (meta data)
% Plot source positions, receiver positions, orientations
SOFAplotGeometry(Obj);

%% Plot energy-time curve (ETC) in the horizontal plane (measurement data)
figure; % for receiver 1  (left ear)
SOFAplotHRTF(Obj,'EtcHorizontal');
figure; % for receiver 2  (right ear)
SOFAplotHRTF(Obj,'EtcHorizontal',2);

%% Plot magnitude spectra in the median plane (measurement data)
figure;
SOFAplotHRTF(Obj,'MagMedian');

%% Download an example audio file (wav)
% Optionally use 'audioread' to load a local wav file
% The following command works for Matlab and Octave
[AudioInput, SamplingRate] = audioread(urlwrite('http://piotr.majdak.com/download/soundlib/raw/gaussiannoise2000ms10msfadeinout.wav', 'MyAudioFile.wav'));

%% Spatialize audio file with hrtf file
% Define audio object, hrtf object, azimuth, elevation
azi=60; ele=0;
[AudioOutput] = SOFAspat(AudioInput, Obj, azi, ele);

%% Play spatialized audio file
% Play; headphones are recommended
sound(AudioOutput, SamplingRate);

%% Spatialize by finding filter, and convolving
% Define azimuth, elevation
azi=60; ele=0;
% Find source position in meta data, get index (= direction)
idx = find(Obj.SourcePosition(:,1)==azi & Obj.SourcePosition(:,2)==ele);
% convolve audio object with hrtf filter with determined index
AudioOutput = [conv(squeeze(Obj.Data.IR(idx, 1, :)), AudioInput) ...
               conv(squeeze(Obj.Data.IR(idx, 2, :)), AudioInput)];
% Play; headphones are recommended
sound(AudioOutput, SamplingRate);

%% Create a new SOFA file and save it
% Get an emtpy object, for convention 'SimpleFreeFieldHRIR'
clear;
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

% Define positions (azimuth, elevation, radius)
azi=-135:45:180;
ele=-30:30:60;
radius=1.2;

% Define dimensions
N=256; % signal length
M=length(azi) * length(ele); % number of measurements

% Create 'data': We create a single impulse after 100 samples; signal length: 256 samples
IR=[zeros(100,1); 1; zeros(N-100-1,1)];

% Prepare size of data in SOFA object
Obj.Data.IR = NaN(M,2,N);

% Fill SOFA object with data and source positions
ii=1;
for aa=1:length(azi)
  for ee=1:length(ele)
    Obj.Data.IR(ii,1,:)=IR;
    Obj.Data.IR(ii,2,:)=IR;
    Obj.SourcePosition(ii,:)=[azi(aa) ele(ee) 1];
    ii=ii+1;
  end
end

% Update dimensions
Obj=SOFAupdateDimensions(Obj);

% Fill some global attributes
Obj.GLOBAL_ListenerShortName = 'Subject 42';
Obj.GLOBAL_History = 'this SOFA file has been created with a demo script';
Obj.GLOBAL_DatabaseName = 'Test Database';
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA Toolbox';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
Obj.GLOBAL_Organization = 'Acoustics Research Institute';
Obj.GLOBAL_AuthorContact = 'michael.mihocic@oeaw.ac.at';
Obj.GLOBAL_Comment = 'This SOFA objects provides first steps to create another object which makes more sense...';

% Save the SOFA file
Obj=SOFAsave('TestSOFAfile.sofa', Obj);

%% SOFA strings
disp(Obj.API.S); % show current (max) strings dimension

% Add strings to SOFA object
Obj = SOFAaddVariable(Obj, 'MySingleString', 'IS', 'This string describes something that cannot described in the other meta data...'); % add string (having dimension 1xS)
myStrings = num2str((1:M)', 'Describing measurement %d'); % one string for each measurement, saved as character array
Obj = SOFAaddVariable(Obj, 'StringForEachMeasurement', 'MS', myStrings); % add string (having dimension MxS)

% Update dimensions
Obj=SOFAupdateDimensions(Obj);
disp(Obj.API.S); % show current (max) strings dimension (updated)

%% Other types than HRTFs: SimpleHeadphoneIR
% Download and load SOFA file
Obj = SOFAload(urlwrite('https://sofacoustics.org/data/examples/SimpleHeadphoneIR_1.0.sofa', 'MyHpIRFile.sofa'));
%% Inspect the data
disp('Inspect the data structure of the SimpleHeadphoneIR example');
keyboard; % Continue with dbcont

% Plot amplitude spectra
% prepare figure to add multiple measurements
figure; hold on;

% Plot all measurements for left ear
for ii=1:Obj.API.M
  plot(20*log10(abs(fft(squeeze(Obj.Data.IR(ii,1,:)),Obj.Data.SamplingRate))) );
  leg{ii}=['Left #' num2str(ii)]; % legend
end

% Plot all measurements for right ear [-20 dB]
for ii=1:Obj.API.M
  plot(20*log10(abs(fft(squeeze(Obj.Data.IR(ii,2,:)),Obj.Data.SamplingRate)))-20 );
  leg{ii+Obj.API.M}=['Right #' num2str(ii)]; % legend
end
xlim([-200 18200]); % range for x axis
legend(leg);
title('Amplitude Spectra of Headphones Measurements: Left, Right [-20 dB]')
xlabel('Frequency (Hz)'); ylabel('Amplitude (dB)');

%% Other types than HRTFs: SingleRoomSRIR
% Download and load SOFA file
Obj = SOFAload('db://database/pan-ar/2_IR_A.sofa');
% Inspect the data
disp('Inspect the data structure of the SingleRoomSRIR example');
keyboard; % Continue with dbcont
% Plot room geometry for room impulse response measurement
SOFAplotGeometry(Obj); % Considers shoebox as the room
view(45,30); % adapt view angle
