function SOFAinfo(Obj)
%SOFAinfo - Display informatio about the SOFA object
%   Usage: SOFAinfo(Obj)
%
%   SOFAinfo(Obj) gathers the mandatory information about the 
%   SOFA object Obj and display it in a narrative form. 
%   For SimpleFreeFieldHRIR more details are displayed.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Global fields

disp('### SOFAinfo: ###')
disp([Obj.GLOBAL_Conventions ' ' Obj.GLOBAL_Version])
disp(['Conventions: ' Obj.GLOBAL_SOFAConventions ' ' Obj.GLOBAL_SOFAConventionsVersion])
disp(['API:         ' Obj.GLOBAL_APIName ' ' Obj.GLOBAL_APIVersion])
disp(['Data Type:   ' Obj.GLOBAL_DataType])
disp(['Room Type:   ' Obj.GLOBAL_RoomType])
disp(['Date Created:   ' Obj.GLOBAL_DateCreated])
disp(['Date Modified:  ' Obj.GLOBAL_DateModified])
disp(['Author Contact: ' Obj.GLOBAL_AuthorContact])
disp(['Organization:   ' Obj.GLOBAL_Organization])
disp(['License:     ' Obj.GLOBAL_License])
disp(['Title:       ' Obj.GLOBAL_Title])


%% Conventions dependant fields

switch Obj.GLOBAL_SOFAConventions

  case 'SimpleFreeFieldHRIR'
disp(['Datebase Name:  ' Obj.GLOBAL_DatabaseName])
      disp(' ')
      disp(['Measurement details (' Obj.GLOBAL_SOFAConventions '):']);
      disp('--------------------');
      disp(['Number of azimuth angles:   ' num2str(length(unique(Obj.SourcePosition(:,1))))]);
      disp(['Number of elevation angles: ' num2str(length(unique(Obj.SourcePosition(:,2))))]);
      disp(['Number of radii:            ' num2str(length(unique(Obj.SourcePosition(:,3))))]);
      if isfield(Obj.Data,'SamplingRate')
          disp(['Sampling Rate: ' num2str(Obj.Data.SamplingRate) ' ' Obj.Data.SamplingRate_Units])
        
      end
      disp(['Listener Short Name: ' Obj.GLOBAL_ListenerShortName])

    if size(Obj.ListenerPosition,1)==1
        disp(['Listener Position: ', num2str(Obj.ListenerPosition) ' ' Obj.ListenerPosition_Units]);
    else
        disp([num2str(size(Obj.ListenerPosition,1)) ' Listener Positions']); 
        disp(['   from ' num2str(Obj.ListenerPosition(1,:)) ' [' Obj.ListenerPosition_Units ']']); 
        disp(['   to   ' num2str(Obj.ListenerPosition(end,:)) ' [' Obj.ListenerPosition_Units ']']);
    end

    if size(Obj.ListenerView,1)==1
        disp(['Listener View:     ', num2str(Obj.ListenerView) ' ' Obj.ListenerView_Units]);
    else
        disp([num2str(size(Obj.ListenerPosition,1)) ' Listener Views']); 
        disp(['   from ' num2str(Obj.ListenerView(1,:)) ' [' Obj.ListenerView_Units ']']); 
        disp(['   to   ' num2str(Obj.ListenerView(end,:)) ' [' Obj.ListenerView_Units ']']);
    end    

    if size(Obj.ReceiverPosition,1)==1
        disp(['Receiver Position: ', num2str(Obj.ReceiverPosition) ' ' Obj.ReceiverPosition_Units]);
    else
        disp([num2str(size(Obj.ReceiverPosition,1)) ' Receiver Positions']); 
        disp(['   from ' num2str(Obj.ReceiverPosition(1,:)) ' [' Obj.ReceiverPosition_Units ']']); 
        disp(['   to   ' num2str(Obj.ReceiverPosition(end,:)) ' [' Obj.ReceiverPosition_Units ']']);
    end
    disp(['Receiver Position Type: ' Obj.ReceiverPosition_Type])

    if size(Obj.SourcePosition,1)==1
        disp(['Source Position:   ', num2str(Obj.SourcePosition) ' ' Obj.SourcePosition_Units]);
    else
        disp([num2str(size(Obj.SourcePosition,1)) ' Source Positions']); 
        disp(['   from ' num2str(Obj.SourcePosition(1,:)) ' [' Obj.SourcePosition_Units ']']); 
        disp(['   to   ' num2str(Obj.SourcePosition(end,:)) ' [' Obj.SourcePosition_Units ']']);
    end
    disp(['Source Position Type: '   Obj.SourcePosition_Type])
  otherwise
    % fill with conventions and disp commands if you feel motivated
end

disp('###    end    ###')
