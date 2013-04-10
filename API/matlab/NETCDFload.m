function [Obj] = NETCDFload(filename,varargin)
%NETCDFLOAD
%   [varName,varContent] = NETCDFload(filename,ReturnType) reads all data (no metadata) from
%   a SOFA file.
%
%   filename specifies the SOFA file from which the data is read.

% SOFA API - function matlab/NETCDFload
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% Global definitions
glob='GLOBAL_';
globid=netcdf.getConstant('GLOBAL');

try
	var='opening file';
	ncid = netcdf.open(filename,'NC_NOWRITE');% open file
	var='inquirying data';
	[numdims,numvars,numglob]  = netcdf.inq(ncid); % get number of anything
		
%% Load global attributes
	for ii=0:numglob-1
    var = netcdf.inqAttName(ncid,globid,ii);
		Obj.(['GLOBAL_' var]) = netcdf.getAtt(ncid,globid,var);
	end
	
%% Load dimensions
	dimids=netcdf.inqDimIDs(ncid);
	dims=cell(numdims,1);
	for ii=0:numdims-1
    [var,len] = netcdf.inqDim(ncid,dimids(ii+1));
		Obj.(var) = len;
		dims{ii+1}=var;
	end
	
%% Load variables and their attributes

	varids=netcdf.inqVarIDs(ncid);
	for ii=0:numvars-1
    [var,~,~,natts] = netcdf.inqVar(ncid,varids(ii+1));	
		if isempty(cell2mat(strfind(dims,var)))	% don't load the data for dimension variables
			data=netcdf.getVar(ncid,varids(ii+1));
			if strfind(var,'Data.'), Obj.Data.(var(6:end))=data; else	Obj.(var)=data; end
		end
		if natts
			for jj=0:natts-1
				att = netcdf.inqAttName(ncid,varids(ii+1),jj);
				attval = netcdf.getAtt(ncid,varids(ii+1),att);
				if strfind(var,'Data.'), Obj.Data.([var(6:end) '_' att])=attval; else Obj.([var '_' att])=attval; end
			end
		end	
	end
	
catch ME
	netcdf.abort(ncid);
	error(['Error processing ' var ' (line ' num2str(ME.stack.line) ')' 10 ...
					'Error message: ' ME.message]);
end

netcdf.close(ncid);
