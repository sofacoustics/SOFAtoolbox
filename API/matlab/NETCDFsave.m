function NETCDFsave(filename,Obj,Var,DataVar,Compression)
%NETCDFSAVE
%   NETCDFsave(filename,Dataset,Compression) saves all data and metadata to
%   a SOFA file.

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
	var='file creation';
	ncid = netcdf.create(filename,'netcdf4');

%% Save global attributes

	for ii=1:length(f)
		if ~isempty(strfind(f{ii},glob))
			var=f{ii};
			id = netcdf.getConstant('GLOBAL');
      netcdf.putAtt(ncid,id,var(strfind(var,glob)+length(glob):end),Obj.(var));
		end
	end
	
%% Define dimensions

	dimid=nan(size(dims));
	dimsize=nan(size(dims));
	for ii=1:length(dims)
		var=dims(ii);
		if isfield(Obj, var)			
			dimid(ii) = netcdf.defDim(ncid,dims(ii),Obj.(var)); 
			dimsize(ii)=Obj.(var);
		end
	end
	netcdf.endDef(ncid);

%% Save dimension variables

	for ii=1:length(dimsize)
		var=dims(ii);
		if ~isnan(dimid(ii))
			VarId = netcdf.defVar(ncid,var,netcdf.getConstant('NC_FLOAT'),dimid(ii));
			netcdf.putVar(ncid,VarId,1:Obj.(var));
			for jj=1:length(f)
				if ~isempty(strfind(f{jj},[var '_']))
					netcdf.putAtt(ncid,VarId,f{jj}(strfind(f{jj},[var '_'])+length([var '_']):end),Obj.(f{jj}));
				end
			end
		else
% 			disp(['DIM skipped: ' dims(ii)]);
		end
	end

%% Save other metadata variables and their attributes

	for ii=1:length(fv)
		var=fv{ii};
		if isempty(strfind(var,'_'))	% skip all attributes	
			ids=cell2mat(regexp(dims,cellstr((Var.(var))')));
			varId = netcdf.defVar(ncid,var,netcdf.getConstant('NC_FLOAT'),dimid(ids));	
			netcdf.putVar(ncid,varId,Obj.(var));
			for jj=1:length(f)
				if ~isempty(strfind(f{jj},[var '_']))
					netcdf.putAtt(ncid,varId,f{jj}(strfind(f{jj},[var '_'])+length([var '_']):end),Obj.(f{jj}));
				end
			end		
		end
	end

%% Save data variables and their attributes
	
	for ii=1:length(fd)
		var=fd{ii};
		if isempty(strfind(var,'_'))	% skip all attributes				
			ids=cell2mat(regexp(dims,cellstr((DataVar.(var))')));
			varId = netcdf.defVar(ncid,['Data.' var],netcdf.getConstant('NC_FLOAT'),dimid(ids));	
			netcdf.putVar(ncid,varId,Obj.Data.(var));
			for jj=1:length(fod)
				if ~isempty(strfind(fod{jj},[var '_']))
					netcdf.putAtt(ncid,varId,fod{jj}(strfind(fod{jj},[var '_'])+length([var '_']):end),Obj.Data.(fod{jj}));
				end
			end		
		end
	end
	
catch ME
	if ~strcmp(ME.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
		netcdf.close(ncid);
	end
	error(['Error processing ' var ' (line ' num2str(ME.stack.line) ')' 10 ...
					'Error message: ' ME.message]);
% 	throw(ME);
end

netcdf.close(ncid);
	
