%demo_GeneralTF - Demonstrates the usage of the GeneralTF conventions.

% #Author: Michael Mihocic
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% 
% SOFA Toolbox - demo script
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Get an empy conventions structure
% clear;
conventions='GeneralTF';
disp(['Creating SOFA file with ' conventions 'conventions...']);
Obj = SOFAgetConventions(conventions);

%% Fill data...
Obj.ReceiverPosition=[0,0.09,0];
Obj.SourcePosition=[1,0,0;0,1,0;-1,0,0;0,-1,0];
% Fill some random data 
Obj.Data.Real=[-0.0055, 0.0003, 0.0001, -0.0005; 0.0009, -0.0014, 0.0010, -0.0024; -0.0055, 0.0003, 0.0001, -0.0005; 0.0009, -0.0014, 0.0010, -0.0024];
Obj.Data.Imag=[0.0009, -0.0014, 0.0010, -0.0024; -0.0055, 0.0003, 0.0001, -0.0005; 0.0009, -0.0014, 0.0010, -0.0024; -0.0055, 0.0003, 0.0001, -0.0005];  

%% Update dimensions
Obj=SOFAupdateDimensions(Obj);

%% Fill with attributes
Obj.GLOBAL_ListenerShortName = 'dummy';
Obj.GLOBAL_History = 'created with a demo script';
Obj.GLOBAL_DatabaseName = 'none';
Obj.GLOBAL_ApplicationName = 'Demo of the SOFA Toolbox';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
Obj.GLOBAL_Organization = 'Acoustics Research Institute';
Obj.GLOBAL_AuthorContact = 'michael.mihocic@oeaw.ac.at';

%% save the SOFA file
SOFAfn=fullfile(SOFAdbPath,'sofatoolbox_test',[conventions '.sofa']);
disp(['Saving:  ' SOFAfn]);
Obj=SOFAsave(SOFAfn, Obj);
