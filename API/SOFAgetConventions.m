function Obj = SOFAgetConventions(sofaconventions)
%SOFAGETVARIABLES
%   Obj = SOFAgetConventions(sofaconventions) returns a list of variables.

% SOFA API - function SOFAgetConventions
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% Last Update: Michael Mihocic, 04.2013


%% ---------------------- return specified variables ----------------------
%% SimpleFreeFieldHRTF
switch sofaconventions
    case 'SimpleFreeFieldHRTF'
       
% Define general information about the file
        Obj.Conventions='SOFA';
        Obj.SOFAConventions=sofaconventions;
        Obj.SOFAConventionVersion='0.0.1';
        Obj.APIName='Matlab ARI';
        Obj.APIVersion='0.0.1';
        Obj.ApplicationName='';
        Obj.ApplicationVersion='';
        Obj.AuthorContact='';
        Obj.Licence='CC BY-SA 3.0';
        Obj.Organization='';
        Obj.DatabaseName='';
        Obj.SubjectID='';
        Obj.R=NaN;
        Obj.E=NaN;
        Obj.N=NaN;
        Obj.M=NaN;
        Obj.C=3;
        Obj.RLongName='number of receivers';
        Obj.ELongName='number of emitters';
        Obj.MLongName='number of measurements';
        Obj.CLongName='coordinate triplet';
        Obj.RoomType='free field';

% Data-type        
        Obj.Data.IR = zeros(3,2,1);
        Obj.DataType = 'IR';
        Obj.Data.SamplingRate = 48000;
        Obj.Data.SamplingRateUnits = 'hertz';
        Obj.NLongName='time';
        Obj.NUnits='samples';

%  Listener
        Obj.ListenerPosition = NaN(1,3);
        Obj.ListenerPositionType = 'cartesian';
        Obj.ListenerPositionUnits = 'meter';
        Obj.ListenerUp = NaN(1,3);
        Obj.ListenerUpType = 'cartesian';
        Obj.ListenerUpUnits = 'meter';
        Obj.ListenerView = NaN(1,3);
        Obj.ListenerViewType = 'cartesian';
        Obj.ListenerViewUnits = 'meter';
        Obj.ListenerRotation = NaN(1,3);
        Obj.ReceiverRotationType = 'din9300';
        Obj.ReceiverRotationUnits = 'degrees';

% Receivers        
%         Obj.ReceiverPosition = [0 -0.09 0; 0 0.09 0];
        Obj.ReceiverPosition = NaN(2,3);
        Obj.ReceiverPositionType = 'cartesian';
        Obj.ReceiverPositionUnits = 'meter';        
        
% Source        
        Obj.SourcePosition = NaN(1,3);
        Obj.SourcePositionType = 'cartesian';
        Obj.SourcePositionUnits = 'meter';          
        Obj.SourceUp = NaN(1,3);
        Obj.SourceUpType = 'cartesian';
        Obj.SourceUpUnits = 'meter'; 
        Obj.SourceView = NaN(1,3);
        Obj.SourceViewType = 'cartesian';
        Obj.SourceViewUnits = 'meter'; 
        
% Emitters
        Obj.EmitterPosition = NaN(1,3);
        Obj.EmitterPositionType = 'cartesian';
        Obj.EmitterPositionUnits = 'meter';    

% Transmitters
        Obj.TransmitterPosition = NaN(1,3);
        Obj.TransmitterPositionType = 'cartesian';
        Obj.TransmitterPositionUnits = 'meter';          
 
    case 'SimpleDRIRMicArray'
%% SimpleDRIRMicArray  

        Obj.Conventions='SOFA';
        Obj.SOFAConventions=sofaconventions;
        Obj.SOFAConventionVersion='0.0.1';
        Obj.APIName='Matlab ARI';
        Obj.APIVersion='0.0.1';
        Obj.ApplicationName='';
        Obj.ApplicationVersion='';
        Obj.AuthorContact='';
        Obj.Licence='CC BY-SA 3.0';
        Obj.Organization='';
        Obj.DatabaseName='';
        Obj.SubjectID='';
        Obj.R=NaN;
        Obj.E=NaN;
        Obj.N=NaN;
        Obj.M=NaN;
        Obj.C=3;
        Obj.RLongName='number of receivers';
        Obj.ELongName='number of emitters';
        Obj.MLongName='number of measurements';
        Obj.CLongName='coordinate triplet';
        Obj.RoomType='dae';

% Data-type        
        Obj.Data.IR = zeros(3,2,1);
        Obj.DataType = 'IR';
        Obj.Data.SamplingRate = 48000;
        Obj.Data.SamplingRateUnits = 'hertz';
        Obj.NLongName='time';
        Obj.NUnits='samples';

%  Listener
        Obj.ListenerPosition = NaN(1,3);
        Obj.ListenerPositionType = 'cartesian';
        Obj.ListenerPositionUnits = 'meter';
        Obj.ListenerUp = NaN(1,3);
        Obj.ListenerUpType = 'cartesian';
        Obj.ListenerUpUnits = 'meter';
        Obj.ListenerView = NaN(1,3);
        Obj.ListenerViewType = 'cartesian';
        Obj.ListenerViewUnits = 'meter';

% Receivers        
%         Obj.ReceiverPosition = [0 -0.09 0; 0 0.09 0];
        Obj.ReceiverPosition = NaN(2,3);
        Obj.ReceiverPositionType = 'spherical';
        Obj.ReceiverPositionUnits = 'degrees';        
        
% Source        
        Obj.SourcePosition = NaN(1,3);
        Obj.SourcePositionType = 'cartesian';
        Obj.SourcePositionUnits = 'meter';          
        Obj.SourceUp = NaN(1,3);
        Obj.SourceUpType = 'cartesian';
        Obj.SourceUpUnits = 'meter'; 
        Obj.SourceView = NaN(1,3);
        Obj.SourceViewType = 'cartesian';
        Obj.SourceViewUnits = 'meter'; 
        
% Emitters
        Obj.EmitterPosition = NaN(1,3);
        Obj.EmitterPositionType = 'cartesian';
        Obj.EmitterPositionUnits = 'meter';    

% Transmitters
        Obj.TransmitterPosition = NaN(1,3);
        Obj.TransmitterPositionType = 'cartesian';
        Obj.TransmitterPositionUnits = 'meter';     
    
        Obj.RoomDAEFileName = 'room';
        Obj.RoomDAEFileNameDescription = 'a room';
    otherwise
        
end

end %of function