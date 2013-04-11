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

%% Check if SOFA conventions
if ~isfield(Obj,'GLOBAL_Conventions'), error('File is not a valid SOFA file'); end
if ~strcmp(Obj.GLOBAL_Conventions,'SOFA'), error('File is not a valid SOFA file'); end
if ~isfield(Obj,'GLOBAL_SOFAConventions'), error('Information about SOFA conventions is missing'); end
try
	X=SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'m');
catch ME
	error('Unsupported SOFA conventions');
end

%% Check if Data present and if correct DataType
if ~isfield(Obj,'Data'), error('Data is missing'); end
if ~isfield(Obj,'GLOBAL_DataType'), error('DataType is missing'); end
if ~strcmp(Obj.GLOBAL_DataType,X.GLOBAL_DataType), error('Wrong DataType'); end
f=fieldnames(X.Data);
for ii=1:length(f)
	if ~isfield(Obj.Data,f{ii}), error(['Data.' f{ii} ' is missing']); end
end

%% Check if mandatory variables are present
f=fieldnames(X);
for ii=1:length(f)
	if ~isfield(Obj,f{ii}), error(['' f{ii} ' is missing']); end
end

Obj=SOFAupdateDimensions(Obj);