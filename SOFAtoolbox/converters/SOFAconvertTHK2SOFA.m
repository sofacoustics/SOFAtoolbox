function Obj = SOFAconvertTHK2SOFA(miroObj)
%SOFAconvertTHK2SOFA - converts from miroObj to SOFA format
%   Usage: OBJ=SOFAconvertTHK2SOFA(miroObj)
% 
%   SOFAconvertTHK2SOFA(miroObj) converts the HRIRs, BRIRs, and DRIRs (VariSphear array measurements) described in miroObj to SOFA. miroObj is the miro object saved at the Technische Hochschule Koeln, provided by Benjamin Bernschuetz.
%
%   Input parameters:
%     miroObj : HRTF data in miro format
% 
%   Output parameters:
%     Obj : New SOFA object (SOFA format)
% 
%   Reference to the source format: http://www.audiogroup.web.th-koeln.de/FILES/miro_documentation.pdf
%   Reference to the source coordinate system: http://www.audiogroup.web.th-koeln.de/SOFiA_wiki/COORDINATES.html

% #Author: Tim Lübeck, TH Köln (2018)
% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%%
%if isoctave,
if exist('OCTAVE_VERSION','builtin') ~= 0
    error(['Octave is not able to convert THK to SOFA, use Matlab instead.']);
end

%% Check if miro.m class file is available, if not download file from server
% miro class might also be included in other toolboxes, eg. AKtools
if exist('miro','class') ~= 8
    % download miro.m
    disp('Downloading miro.m from TH Köln server...')
    url = 'http://audiogroup.web.th-koeln.de/FILES/miro.m'; 
    basepath=which('SOFAstart');
    basepath=basepath(1:end-12); % Kill the function name from the path.
    target=[basepath filesep 'helpers' filesep 'miro.m'];
    websave (target,url);
end

%% Get an empty conventions structure depending on miro format
if ( strcmp (miroObj.type,'HRIR') )
    Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
    Obj.GLOBAL_Title = 'HRIR';
    Obj.GLOBAL_RoomType = 'free field';
    typeFlag = 0; 
elseif ( strcmp (miroObj.type,'BRIR') )
    Obj = SOFAgetConventions('MultiSpeakerBRIR');
    Obj.GLOBAL_Title = 'BRIR';
    Obj.GLOBAL_RoomType = 'reverberant';
    typeFlag = 1; 
else %( strcmp (miroObj.type,'MICARRAY') )
    Obj = SOFAgetConventions('SingleRoomDRIR');
    Obj.GLOBAL_Title = 'DRIR';
    Obj.GLOBAL_RoomType = 'reverberant';
    typeFlag = 2; 
end

%% Fill global attributes
Obj.GLOBAL_ListenerShortName = miroObj.name;
Obj.GLOBAL_AuthorContact = miroObj.contact;
Obj.GLOBAL_Comment = [miroObj.context,' / Sampling Grid: ',miroObj.quadGrid];
Obj.GLOBAL_History = SOFAappendText(Obj,'GLOBAL_History','Converted from the miro file format');
Obj.GLOBAL_License = 'CC 3.0 BY-SA';
Obj.GLOBAL_Organization = 'Technische Hochschule Koeln, Germany';
Obj.GLOBAL_Author = miroObj.engineer;
Obj.GLOBAL_Origin = 'http://audiogroup.web.th-koeln.de';
Obj.GLOBAL_DateCreated = datestr(datenum(miroObj.date),'yyyy-mm-dd HH:MM:SS');
Obj.GLOBAL_DatabaseName='THK';
Obj.GLOBAL_ListenerDescription = miroObj.microphone;
Obj.GLOBAL_ReceiverDescription = [miroObj.microphone '; ' miroObj.micPreamp];
Obj.GLOBAL_SourceDescription = miroObj.source;
Obj.GLOBAL_EmitterDescription = miroObj.source;
Obj.GLOBAL_RoomDescription = [miroObj.location,' / avgAirTemp: ',num2str(miroObj.avgAirTemp),' / avgRelHumidity: ',num2str(miroObj.avgRelHumidity)];

%% Set miroObj to degree mode
miroObj.shutUp = 1;
miroObj = setDEG(miroObj);

%% Get miroObject data and convert to SOFA structure
%M: number of measurements; 
%R: number of receivers; 
%N: number of data samples describing one measurement. Data is a function of N;
%E: number of emitters;  
%C: coordinate dimension, always three with the meaning
%FIR : [M R N];

if (typeFlag == 0 || typeFlag == 1)  %BRIR or HRIR
    irChOne = zeros(miroObj.returnTaps,miroObj.nIr);
    irChTwo = zeros(miroObj.returnTaps,miroObj.nIr);
    for channel = 1 : miroObj.nIr
        IR = getIR( miroObj, channel );
        irChOne(:, channel) = IR(:,1);  %[N M]
        irChTwo(:, channel) = IR(:,2);  %[N M]   
    end 
    if (typeFlag == 0) %HRIR
        Obj.Data.IR = irChOne; % irChOne is [N M]
        Obj.Data.IR(:,:,2) = irChTwo;
        Obj.Data.IR = shiftdim(Obj.Data.IR,1); % convert from [N M R] to [M R N]
    else % BRIR
        Obj.Data.IR = zeros(size(irChOne,2), 2, 1, size( irChTwo,1)); %[N 2 1 M]
        Obj.Data.IR(:,1,:) = shiftdim(shiftdim(irChOne,-2),3); 
        Obj.Data.IR(:,2,:) = shiftdim(shiftdim(irChTwo,-2),3);
    end
else % DRIR
    irData = zeros(miroObj.returnTaps,miroObj.nIr);
    for channel = 1 : miroObj.nIr
        irData(:,channel) = miroObj.getIR(channel);
    end
    Obj.Data.IR             = zeros(1, miroObj.nIr ,miroObj.returnTaps); %[M R N]
    for R = 1 : miroObj.nIr
        Obj.Data.IR(1,R,:) = irData(:,R);
    end
    Obj.Data.Delay = zeros(1,miroObj.nIr);
end
Obj.Data.SamplingRate = miroObj.fs;

%% Fill the mandatory source emitter variables
Obj.ListenerPosition = [0 0 0]; % for BRIR and HRIR listener in center
Obj.ReceiverPosition = [0 +miroObj.radius 0; 0 -miroObj.radius 0];  % for HRIR and BRIR ears as receiver 

if (typeFlag == 0)%HRIR
    Obj.ListenerView = [1 0 0];
    Obj.ListenerUp = [0 0 1];
    Obj.SourcePosition = [...
        miroObj.azimuth' ... % azimuth angle in a range of (0-360°(. Whereas AZ=0° is defined to be the front direction and AZ=180° to be the rear direction.
        90-miroObj.elevation' ... % elevation angle in range of (0-180°). EL=0 points upwards, EL=90° points to the horizontal plane and EL=180° points downwards.
        miroObj.sourceDistance*ones(size(miroObj.azimuth'))]; % radius in meters
elseif (typeFlag == 1) %BRIR
    Obj.SourcePosition = [0 0 0]; % default edit manually!
    Obj.EmitterPosition = [miroObj.sourceDistance 0 0];   % default position is center, otherwise define manually
    Obj.EmitterPosition_Type  = 'cartesian';
    Obj.EmitterUp = [0 0 1];
    Obj.EmitterView = [-1 0 0];
    Obj.ListenerView = [miroObj.azimuth',  ...          % see HRIR definitions
                        90-miroObj.elevation',  ...
                        zeros(size(miroObj.azimuth'))]; %miroObj.sourceDistance*ones(size(miroObj.azimuth'))]; 
    Obj.ListenerView_Type = 'spherical';
    Obj.ListenerView_Units = 'degree, degree, metre';
    Obj.ListenerUp = [0 0 1]; 
else %DRIR
    Obj.SourcePosition  = [1,0,0]; % default edit manually!
    Obj.EmitterPosition = [0,0,0];
    Obj.ListenerPosition    = [0 0 0];
    Obj.ListenerView        = [1 0 0];
    Obj.ListenerUp          = [0 0 1];
    Obj.ListenerView_Type = 'spherical';
    Obj.ListenerView_Units = 'degree, degree, metre';
    Obj.ReceiverPosition = [...
                            miroObj.azimuth' ... 
                            90-miroObj.elevation' ... 
                            miroObj.radius*ones(size(miroObj.azimuth'))]; 
    Obj.ReceiverPosition_Type = 'spherical';
    Obj.ReceiverPosition_Units = 'degree, degree, metre';
end

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);
