function NETCDFsave(filename,Obj,Var,DataVar,Compression)
%NETCDFSAVE
%   NETCDFsave(filename,Dataset,Compression) saves all data and metadata to
%   a SOFA file.
% dim: [M,R,N,E,C,Q]

% SOFA API - function matlab/NETCDFsave
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% Piotr Majdak, 9.4.2013



%% Global definitions
glob='GLOBAL_';
dims='MRNECQ';
f=fieldnames(Obj);
fv=fieldnames(Var);
fd=fieldnames(DataVar);
fod=fieldnames(Obj.Data);

try 
	ncid = netcdf.create(filename,'netcdf4');

%% Save global attributes

	for ii=1:length(f)
		if ~isempty(strfind(f{ii},glob))
			id = netcdf.getConstant('GLOBAL');
      netcdf.putAtt(ncid,id,f{ii}(strfind(f{ii},glob)+length(glob):end),Obj.(f{ii}));
		end
	end
	
%% Define dimensions

	dimid=nan(size(dims));
	dimsize=nan(size(dims));
	for ii=1:length(dims)
		if isfield(Obj, dims(ii))
			dimid(ii) = netcdf.defDim(ncid,dims(ii),Obj.(dims(ii))); 
			dimsize(ii)=Obj.(dims(ii));
		end
	end
	netcdf.endDef(ncid);

%% Save dimension variables

	for ii=1:length(dimsize)
		if ~isnan(dimid(ii))
			disp(['DIM: ' dims(ii)]);
			VarId = netcdf.defVar(ncid,dims(ii),netcdf.getConstant('NC_FLOAT'),dimid(ii));
			netcdf.putVar(ncid,VarId,1:Obj.(dims(ii)));
			for jj=1:length(f)
				if ~isempty(strfind(f{jj},[dims(ii) '_']))
					netcdf.putAtt(ncid,VarId,f{jj}(strfind(f{jj},[dims(ii) '_'])+length([dims(ii) '_']):end),Obj.(f{jj}));
				end
			end
		else
			disp(['DIM skipped: ' dims(ii)]);
		end
	end

%% Save other metadata variables and their attributes

	for ii=1:length(fv)
		if isempty(strfind(fv{ii},'_'))	% skip all attributes	
			disp(['VAR: ' fv{ii}]);
			ids=cell2mat(regexp(dims,cellstr((Var.(fv{ii}))')));
			varId = netcdf.defVar(ncid,fv{ii},netcdf.getConstant('NC_FLOAT'),dimid(ids));	
			netcdf.putVar(ncid,varId,Obj.(fv{ii}));
			for jj=1:length(f)
				if ~isempty(strfind(f{jj},[fv{ii} '_']))
					netcdf.putAtt(ncid,varId,f{jj}(strfind(f{jj},[fv{ii} '_'])+length([fv{ii} '_']):end),Obj.(f{jj}));
				end
			end		
		end
	end

%% Save data variables and their attributes
	
	for ii=1:length(fd)
		if isempty(strfind(fd{ii},'_'))	% skip all attributes	
			disp(['DATA: ' fd{ii}]);
			ids=cell2mat(regexp(dims,cellstr((DataVar.(fd{ii}))')));
			varId = netcdf.defVar(ncid,['Data.' fd{ii}],netcdf.getConstant('NC_FLOAT'),dimid(ids));	
			netcdf.putVar(ncid,varId,Obj.Data.(fd{ii}));
			for jj=1:length(fod)
				if ~isempty(strfind(fod{jj},[fd{ii} '_']))
					netcdf.putAtt(ncid,varId,fod{jj}(strfind(fod{jj},[fd{ii} '_'])+length([fd{ii} '_']):end),Obj.Data.(fod{jj}));
				end
			end		
		end
	end
	
catch ME
	if ~strcmp(ME.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
		netcdf.close(ncid);
	end
	error([ME.message ' Error in line ' num2str(ME.stack.line)]);
end

netcdf.close(ncid);
	
