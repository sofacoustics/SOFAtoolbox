function [Obj,Dims] = NETCDFload(filename,flags,varargin)
%% NETCDFLOAD
%   Obj = NETCDFload(filename,'all') reads the SOFA object OBJ with all data from
%   a SOFA file.
%
%   Obj = NETCDFload(filename,'nodata') ignores the Data. variables while
%   reading.
%
%   Obj = NETCDFload(filename,[START COUNT],partialDim) reads only COUNT number of data
%   in dimension partialDim (given as string) beginning with the index START. If START
%   and COUNT are column vectors, then partialDim has to be a string containing the
%   dimension for every column entry.
%
%   [Obj,Dims] = NETCDFload(...) returns the dimension variables found in
%   the file as a string.

% 10.8.2014: string array support. Works for 1D and 2D strings only.
% 
% SOFA API - function netcdf/NETCDFload
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% If we are running octave we have to import the NETCDF namespace, in order to
% run functions like netcdf.getConstant
if exist('OCTAVE_VERSION','builtin')
    import_netcdf;
end

%% Global definitions
glob = 'GLOBAL_';
globid = netcdf.getConstant('GLOBAL');

%% --------------------------- N E T C D F load --------------------------
% Open the NETCDF file
try
	var = 'opening file';
	ncid = netcdf.open(filename,'NC_NOWRITE');      % open file
	var = 'inquirying data';
	[numdims,numvars,numglob]  = netcdf.inq(ncid);  % get number of anything
		
    % ----- GLOBAL ATTRIBUTES --------------------------------------------
	for ii=0:numglob-1
        var = netcdf.inqAttName(ncid,globid,ii);
		Obj.([glob var]) = netcdf.getAtt(ncid,globid,var);
	end
	
    % ----- DIMENSIONS ---------------------------------------------------
	dimids = netcdf.inqDimIDs(ncid);
	dims = cell(numdims,1);    % cell array with dimension names
	startp = zeros(numdims,1); % vector with start of a dimension
	countp = zeros(numdims,1); % vector with the element count in a dimension
	for ii=0:numdims-1
        [var,len] = netcdf.inqDim(ncid,dimids(ii+1));
		Obj.API.(var) = len;
		dims{ii+1} = var;
		startp(ii+1) = 0;
		countp(ii+1) = len;
	end
	Dims = cell2mat(dims)';
	
    % Check the requested measurements
    if isnumeric(flags)
        partialDimRange = flags;
        partialDim = varargin{1};
        for ii=1:length(partialDim)
            if Obj.API.(partialDim(ii))<(sum(partialDimRange(ii,:))-1)
                error('Requested indices exceed measurement count');
            end;
            startp(strfind(Dims,partialDim(ii))) = partialDimRange(ii,1)-1;
            countp(strfind(Dims,partialDim(ii))) = partialDimRange(ii,2);
        end
    end
	

    % ----- VARIABLES + ATTRIBUTES ---------------------------------------
	varids = netcdf.inqVarIDs(ncid);
	for ii=0:numvars-1
        [var,~,vardimids,natts] = netcdf.inqVar(ncid,varids(ii+1));	
        % Data
        if strfind(var,'Data.'),
            if ~strcmp(flags,'nodata')
                dim=fliplr(cell2mat(dims(vardimids+1))');
                Obj.API.Dimensions.Data.(var(6:end))=dim;
                if strfind(dim,'S') % strings
                    if length(dim)>2  % 2D string arrays. ToDo: MdD string arrays
                        data=netcdf.getVar(ncid,varids(ii+1),startp(vardimids+1),countp(vardimids+1));              
                        s=size(data);
                        Obj.Data.(var(6:end))=cell(s(end:-1:2));
                        data=reshape(reshape(data,1,[]),[s(2:end) s(1)]);
                        for jj=1:s(2)
                            for kk=1:s(3)
                                Obj.Data.(var(6:end))(kk,jj)=cellstr(squeeze(data(jj,kk,:))');
                            end
                        end
                    else % 1D string array
                        data=netcdf.getVar(ncid,varids(ii+1),startp(vardimids+1),countp(vardimids+1));
                        Obj.Data.(var(6:end))=cellstr(reshape(reshape(data,1,[]),size(data,2),[]));
                    end
                elseif length(dim)>1
                    Obj.Data.(var(6:end))=permute(netcdf.getVar(ncid,varids(ii+1),startp(vardimids+1),countp(vardimids+1)), length(dim):-1:1); 
                else
                    Obj.Data.(var(6:end))=netcdf.getVar(ncid,varids(ii+1),startp(vardimids+1),countp(vardimids+1));
                end
            end
        % Variables
        else
            data=netcdf.getVar(ncid,varids(ii+1),startp(vardimids+1),countp(vardimids+1));
            dim=fliplr(cell2mat(dims(vardimids+1))');
            Obj.API.Dimensions.(var)=dim;
            if strfind(dim,'S')
                if length(dim)>2  % 2D string arrays. ToDo: MdD string arrays
                    data=netcdf.getVar(ncid,varids(ii+1),startp(vardimids+1),countp(vardimids+1));              
                    s=size(data);
                    Obj.(var)=cell(s(end:-1:2));
                    data=reshape(reshape(data,1,[]),[s(2:end) s(1)]);
                    for jj=1:s(2)
                        for kk=1:s(3)
                            Obj.(var)(kk,jj)=cellstr(squeeze(data(jj,kk,:))');
                        end
                    end
                else % 1D string array
                    data=netcdf.getVar(ncid,varids(ii+1),startp(vardimids+1),countp(vardimids+1));
                    Obj.(var)=cellstr(reshape(reshape(data,1,[]),size(data,2),[]));
                end
            elseif length(dim)>1
                Obj.(var)=permute(data, length(dim):-1:1); 
            else
                Obj.(var)=data;
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
