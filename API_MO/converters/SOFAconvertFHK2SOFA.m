function Obj=SOFAconvertFHK2SOFA(miroObj)
% OBJ=SOFAconvertFHK2SOFA(miroObj) converts the HRTFs described in miroObj
% to SOFA. miroObj is the miro object saved at the Fach-Hochschule Köln
% provided by Benjamin Bernschütz.
% Reference to the source format: http://www.audiogroup.web.fh-koeln.de/ku100hrir.html
% Reference to the source coordinate system: [1] http://code.google.com/p/sofia-toolbox/wiki/COORDINATES

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

if isoctave,
    error(['Octave is not able to convert FHK to SOFA, use Matlab instead.']);
end

%% Get an empty conventions structure
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

%% Fill data with data
Obj.Data.IR = miroObj.irChOne; % irChOne is [N M]
Obj.Data.IR(:,:,2) = miroObj.irChTwo;
Obj.Data.IR = shiftdim(Obj.Data.IR,1); % convert from [N M R] to [M R N]
Obj.Data.SamplingRate = miroObj.fs;

%% Fill with attributes
Obj.GLOBAL_ListenerShortName = miroObj.name;
% Obj.GLOBAL_APIName
% Obj.GLOBAL_ApplicationName
% Obj.GLOBAL_ApplicationVersion
Obj.GLOBAL_AuthorContact = miroObj.contact;
Obj.GLOBAL_Comment = miroObj.context;
Obj.GLOBAL_History = SOFAappendText(Obj,'GLOBAL_History','Converted from the miro file format');
Obj.GLOBAL_License = 'CC 3.0 BY-SA';
Obj.GLOBAL_Organization = 'Fachhochschule Köln, Germany';
Obj.GLOBAL_Author = miroObj.engineer;
Obj.GLOBAL_References = 'Bernschütz, B. (2013). "A Spherical Far Field HRIR/HRTF Compilation of the Neumann KU 100", proceedings of the AIA/DAGA, Meran, Italy';
Obj.GLOBAL_RoomType = 'free field';
Obj.GLOBAL_Origin = 'http://www.audiogroup.web.fh-koeln.de/ku100hrir.html';
Obj.GLOBAL_DateCreated = datestr(datenum(miroObj.date),'yyyy-mm-dd HH:MM:SS');
Obj.GLOBAL_DatabaseName='FHK';
Obj.GLOBAL_Title = 'HRTF';

%% Fill the mandatory variables
Obj.ReceiverPosition = [0 -miroObj.radius 0; 0 +miroObj.radius 0];  
Obj.ListenerPosition = [0 0 0];
Obj.ListenerView = [1 0 0];
Obj.ListenerUp = [0 0 1];
  % From [1]: All angles are in RAD.
  %   AZ denotes the azimutal angle in a range of (0-2pi(. Whereas AZ=0 is defined to be the front direction and AZ=pi to be the rear direction.
  %   EL denotes the elevation angle in a range of (0-pi). EL=0 points upwards, EL=pi/2 points to the horizontal plane and EL=pi points downwards.
  %   r is the radius in meters (if needed/specified). 
Obj.SourcePosition = [...
            rad2deg(miroObj.azimuth') ... % miroObj.azimuth is AZ
            90-rad2deg(miroObj.elevation') ... % miroObj.elevation is EL 
            miroObj.sourceDistance*ones(size(miroObj.azimuth'))]; % miro.sourceDistance is r

Obj.GLOBAL_ListenerDescription = miroObj.microphone;
Obj.GLOBAL_ReceiverDescription = [miroObj.microphone '; ' miroObj.micPreamp];
Obj.GLOBAL_SourceDescription = miroObj.source;
%Obj.GLOBAL_EmitterDescription ='';
Obj.GLOBAL_RoomDescription = miroObj.location;

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);
