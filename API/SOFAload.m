function Obj = SOFAload(filename,flags)
%SOFALOAD 
%   Obj = SOFAload(filename) reads the SOFA object OBJ with all data from
%   a SOFA file.
%
%   Obj = SOFAload(filename,'nodata') ignores the Data. variables and
%   results in the metadata only (variables and attributes).
%
%   Obj = SOFAload(filename,[START COUNT]) reads only COUNT number of 
%		measurements beginning with the index START.
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
if ~exist('flags','var'), flags='all'; end;

%% check file name
filename=SOFAcheckFilename(filename);

%% Load the object
[Obj,Dims]=NETCDFload(filename,flags);

%% Check if SOFA conventions
if ~isfield(Obj,'GLOBAL_Conventions'), error('File is not a valid SOFA file'); end
if ~strcmp(Obj.GLOBAL_Conventions,'SOFA'), error('File is not a valid SOFA file'); end
if ~isfield(Obj,'GLOBAL_SOFAConventions'), error('Information about SOFA conventions is missing'); end
try
	X=SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'m');
catch ME
	error('Unsupported SOFA conventions');
end

%% Check if DataType is present and correct
if ~isfield(Obj,'GLOBAL_DataType'), error('DataType is missing'); end
if ~strcmp(Obj.GLOBAL_DataType,X.GLOBAL_DataType), error('Wrong DataType'); end

%% If data loaded, check if correct data format
if ~strcmp(flags,'nodata'),
	if ~isfield(Obj,'Data'), error('Data is missing'); end
	f=fieldnames(X.Data);
	for ii=1:length(f)
		if ~isfield(Obj.Data,f{ii}), error(['Data.' f{ii} ' is missing']); end
	end
	Obj=SOFAupdateDimensions(Obj);
end

%% Check if mandatory variables are present
f=fieldnames(X);
for ii=1:length(f)
	if ~isfield(Obj,f{ii}), error(['' f{ii} ' is missing']); end
end
