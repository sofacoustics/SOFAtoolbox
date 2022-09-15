function Obj = SOFAconvertSYMARE2SOFA(SYMAREpath,HRIRname)
%SOFAconvertSYMARE2SOFA - converts from SYMARE to SOFA format
%   Usage: Obj = SOFAconvertSYMARE2SOFA(SYMAREpath,HRIRname)
% 
%   SOFAconvertSYMARE2SOFA(SYMAREpath,HRIRname) converts objects from SYMARE database to SOFA format.
%
%   Input parameters:
%     SYMAREpath : Path of the SYMARE directory
%     HRIRname   : HRIR in the <SYMAREpath>/HRIRs/Acoustic directory to be converted
% 
%   Output parameters:
%     Obj : New SOFA object (SOFA format)

% #Author: Piotr Majdak
% #Author: Michael Mihocic: license added, header documentation updated (28.10.2021)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Get an empty conventions structure

Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
Obj.GLOBAL_Title = 'HRIR';
Obj.GLOBAL_RoomType = 'free field';

%% Fill global attributes
if isempty(strfind(HRIRname,'.mat'))
    Obj.GLOBAL_ListenerShortName = strcat('Subj_',HRIRname(end-1:end));
else
    Obj.GLOBAL_ListenerShortName = strcat('Subj_',HRIRname(end-5:end-4));
end

Obj.GLOBAL_AuthorContact ='https://www.morphoacoustics.org/resources.html';
Obj.GLOBAL_Comment = '';
Obj.GLOBAL_History = SOFAappendText(Obj,'GLOBAL_History', ... 
                                    'Converted from the SYMARE database');
Obj.GLOBAL_License = strcat('Creative Commons Attribution-Non', ...
                 'Commercial-ShareAlike 4.0 International Public License');
Obj.GLOBAL_Organization = 'University of Sydney';
Obj.GLOBAL_Author = 'Craig Jin, Anthony Tew, et al.';
Obj.GLOBAL_Origin = 'https://www.morphoacoustics.org/resources.html';
Obj.GLOBAL_DateCreated = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS'));
Obj.GLOBAL_DatabaseName='SYMARE';

%% Get SYMARE data and convert to SOFA structure
%m: number of measurements; 
%R: number of receivers; 
%n: number of data samples describing one measurement. Data is a function of N;
%E: number of emitters;  
%C: coordinate dimension, always three with the meaning
%HRIR Dimensions : mRn;

% Get source positions
load(fullfile(SYMAREpath,'HRIRs','Parameters','azim'));
load(fullfile(SYMAREpath,'HRIRs','Parameters','elev'));
load(fullfile(SYMAREpath,'HRIRs','Parameters','r'));
if ~(length(r) == length(azim))
    r(1:size(azim)) = r;
    r = r';
end
Obj.SourcePosition = [azim*180/pi,elev*180/pi,r];

% Get sampling frequency
load(fullfile(SYMAREpath,'HRIRs','Parameters','fs'));
Obj.Data.SamplingRate = fs;

% Get IRs
load(fullfile(SYMAREpath,'HRIRs','Acoustic',HRIRname)); % gets you hR, hL
HRIR(:,:,1) = hL'; % nm -> mnR
HRIR(:,:,2) = hR'; % nm -> mnR
Obj.Data.IR = permute(HRIR,[1 3 2]); % mnR -> mRn

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

end