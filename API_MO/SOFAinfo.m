function SOFAinfo(Obj)
% SOFAplotGeometry(Obj) plots the geometry found in the Obj.
% 
% SOFAplotGeometry(Obj, index) plots the geometry for the measurements
% given in the index. 

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


% Gather information about the SOFA file and display them



switch Obj.GLOBAL_SOFAConventions
%%
  case 'SimpleFreeFieldHRIR'
    fprintf('\n');
    fprintf('%s\n',Obj.GLOBAL_Title);
    fprintf('%s\n',repmat('=',1,length(Obj.GLOBAL_Title)));
    fprintf('\n');
    fprintf('Anechoic HRTF mesurement done by %s.\n',Obj.GLOBAL_Organization);
    fprintf('\n');
    fprintf('Contact: %s\n',Obj.GLOBAL_AuthorContact);
    fprintf('License: %s\n',Obj.GLOBAL_License);
    fprintf('URL: %s\n',Obj.GLOBAL_Origin);
    fprintf('\n');
    fprintf('Measurement details:\n');
    fprintf('--------------------\n');
    fprintf('Number of azimuth angles:   % 5.0f\n',length(unique(Obj.SourcePosition(:,1))));
    fprintf('Number of elevation angles: % 5.0f\n',length(unique(Obj.SourcePosition(:,2))));
    fprintf('Number of radii:            % 5.0f\n',length(unique(Obj.SourcePosition(:,3))));
    fprintf('Sampling Rate: %.0f %s\n',Obj.Data.SamplingRate,Obj.Data.SamplingRate_Units);
    fprintf('Dummy head: %s\n',Obj.GLOBAL_ListenerShortName);
    fprintf('Loudspeaker: %s\n',Obj.GLOBAL_SourceDescription);
    if size(Obj.ListenerPosition,1)==1
        fprintf('Listener at (%.1f,%.1f,%.1f) %s\n', ...
            Obj.ListenerPosition,Obj.ListenerPosition_Units);
    else
        fprintf(['%.0f listener positions from (%.1f,%.1f,%.1f) %s ', ...
            'to (%.1f,%.1f,%.1f) %s\n'],size(Obj.ListenerPosition,1), ...
            Obj.ListenerPosition(1,:),Obj.ListenerPosition_Units, ...
            Obj.ListenerPosition(end,:),Obj.ListenerPosition_Units);
    end
    if size(Obj.SourcePosition,1)==1
        fprintf('Source at (%.1f,%.1f,%.1f) %s\n', ...
            Obj.SourcePosition,Obj.SourcePosition_Units);
    else
        fprintf(['%.0f source positions from (%.1f,%.1f,%.1f) %s ', ...
            'to (%.1f,%.1f,%.1f) %s\n'],size(Obj.SourcePosition,1), ...
            Obj.SourcePosition(1,:),Obj.SourcePosition_Units, ...
            Obj.SourcePosition(end,:),Obj.SourcePosition_Units);
    end
    fprintf('\n');
   
%%    
  case 'SingleRoomDRIR'    
  otherwise
    error('SOFAConventions not supported');
end
