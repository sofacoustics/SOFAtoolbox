function [Obj,Dims] = NETCDFload(filename,flags)
%NETCDFLOAD
%   Obj = NETCDFload(filename,'all') reads the SOFA object OBJ with all data from
%   a SOFA file.
%
%   Obj = NETCDFload(filename,'nodata') ignores the Data. variables while
%   reading.
%
%   Obj = NETCDFload(filename,[START COUNT]) reads only COUNT number of 
%		measurements beginning with the index START.
%
%   [Obj,Dims] = NETCDFload(...) returns the dimension variables found in
%   the file as a string.


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

%% Open the NETCDF file
try
	var='opening file';
	ncid = netcdf.open(filename,'NC_NOWRITE');% open file
	var='inquirying data';
	[numdims,numvars,numglob]  = netcdf.inq(ncid); % get number of anything
		
%% Load global attributes
	for ii=0:numglob-1
    var = netcdf.inqAttName(ncid,globid,ii);
		Obj.([glob var]) = netcdf.getAtt(ncid,globid,var);
	end
	
%% Load dimensions
	dimids=netcdf.inqDimIDs(ncid);
	dims=cell(numdims,1); % cell array with dimension names
	startp=zeros(numdims,1); % vector with start of a dimension
	countp=zeros(numdims,1); % vector with the element count in a dimension
	for ii=0:numdims-1
    [var,len] = netcdf.inqDim(ncid,dimids(ii+1));
		Obj.(var) = len;
		dims{ii+1}=var;
		startp(ii+1)=0;
		countp(ii+1)=len;
	end
	Dims=cell2mat(dims)';
	
%% Check the requested measurements
if isnumeric(flags)
	if Obj.M<flags(2), error('Requested end index exceeds the measurement count'); end;
	startp(strfind(Dims,'M'))=flags(1)-1;
	countp(strfind(Dims,'M'))=flags(2);
end
	
%% Load variables and their attributes

	varids=netcdf.inqVarIDs(ncid);
	for ii=0:numvars-1
    [var,~,vardimids,natts] = netcdf.inqVar(ncid,varids(ii+1));	
		if isempty(cell2mat(strfind(dims,var)))	% don't load the data for dimension variables			
			if strfind(var,'Data.'),
				if ~strcmp(flags,'nodata')
					data=netcdf.getVar(ncid,varids(ii+1),startp(vardimids+1),countp(vardimids+1));
					Obj.Data.(var(6:end))=data; 
					Obj.Dimensions.Data.(var(6:end))=cell2mat(dims(vardimids+1))';
				end
			else
				data=netcdf.getVar(ncid,varids(ii+1),startp(vardimids+1),countp(vardimids+1));
				Obj.(var)=data; 
				Obj.Dimensions.(var)=cell2mat(dims(vardimids+1))';
			end
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
	if exist('ncid','var'); netcdf.abort(ncid); end;
	for ii=1:length(ME.stack)
		disp(ME.stack(ii));
	end
	error(['Error processing ' var 10 ...
					'Error message: ' ME.message 10 'See also the error stack before']);
end

netcdf.close(ncid);
