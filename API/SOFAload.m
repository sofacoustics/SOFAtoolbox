function Obj = SOFAload(filename,varargin)
%SOFALOAD 
%   Obj = SOFAload(filename,ReturnType) reads all data from a SOFA file.
%   filename specifies the SOFA file from which the data is read.
%   
% SOFA API - function SOFAload
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% Global definitions
dims={'i';'r';'e';'n';'m';'c';'q'}; % dimensions

%% check file name
filename=SOFAcheckFilename(filename);

%% Load the object
[Obj,Dims]=NETCDFload(filename);

% ReturnType = 'struct'; % set default value for ReturnType
% if size(varargin,2)==1
% 	varargin = cellstr(varargin);
% 	ReturnType = varargin{1};
% end
% if ~ischar(ReturnType)
% 	error('ReturnType must be a string.');
% end
% switch ReturnType
%     case 'struct'
%         results = struct; % initialize struct variable
%     case 'cell'
%         
%     otherwise
%         error('ReturnType must be either ''struct'' or ''cell''.');
% end
% 
% %% --------------------------- N E T C D F load ---------------------------
% [varName,varContent]=NETCDFload(filename,'meta');
% for ii=1:length(varName)
%     if strcmp(ReturnType,'struct')
%         results.(varName{ii})=varContent{ii};
%     elseif strcmp(ReturnType,'cell')
%         result{1}=varName{ii};
%         result{2}=varContent{ii};
%         results{ii}=result;
%     end
% end
% 
% [varName,varContent]=NETCDFload(filename,'data');
% for ii=1:length(varName)
%     if strcmp(ReturnType,'struct')
%         results.Data.(varName{ii}(6:end))=varContent{ii};
%     elseif strcmp(ReturnType,'cell')
%         result{1}=['Data' varName{ii}(6:end)];
%         result{2}=varContent{ii};
%         results{length(results)+1}=result;
%     end
% end
% 
