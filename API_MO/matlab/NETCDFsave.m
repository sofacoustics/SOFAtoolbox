function NETCDFsave(filename,Obj,Compression)
%NETCDFSAVE
%   NETCDFsave(filename,Dataset,Compression) saves all data and metadata to
%   a SOFA file.

% SOFA API - function matlab/NETCDFsave
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% Piotr Majdak, 9.4.2013

%% Global definitions
glob='GLOBAL_';
Dims='IRENMCQ'; % dimensions

globid=netcdf.getConstant('GLOBAL');

try 
	var='file creation';
	ncid = netcdf.create(filename,'netcdf4');

%% Save global attributes
f=fieldnames(Obj);

	for ii=1:length(f)
		if ~isempty(strfind(f{ii},glob))
			var=f{ii};
      netcdf.putAtt(ncid,globid,var(strfind(var,glob)+length(glob):end),Obj.(var));
		end
	end
	
%% Define dimensions

	dimid=nan(size(Dims));
	dimsize=nan(size(Dims));
	for ii=1:length(Dims)
		var=Dims(ii);
		if isfield(Obj, var)			
			dimid(ii) = netcdf.defDim(ncid,Dims(ii),Obj.(var)); 
			dimsize(ii)=Obj.(var);
		end
	end
	netcdf.endDef(ncid);

%% Save dimension variables

	for ii=1:length(dimsize)
		var=Dims(ii);
		if ~isnan(dimid(ii))
			VarId = netcdf.defVar(ncid,var,netcdf.getConstant('NC_INT'),dimid(ii));
			netcdf.putVar(ncid,VarId,1:Obj.(var));
			for jj=1:length(f)
				if ~isempty(strfind(f{jj},[var '_']))
					netcdf.putAtt(ncid,VarId,f{jj}(strfind(f{jj},[var '_'])+length([var '_']):end),Obj.(f{jj}));
				end
			end
		end
	end

%% Save other metadata variables and their attributes
Dimensions=rmfield(Obj.Dimensions,'Data');
fv=fieldnames(Dimensions);

	for ii=1:length(fv)
		var=fv{ii};
		if isempty(strfind(var,'_')) && ~strcmp(var,'Dimensions')	% skip all attributes	& dimensions
			ids=cell2mat(regexp(Dims,cellstr((Dimensions.(var))')));
			varId = netcdf.defVar(ncid,var,netcdf.getConstant('NC_DOUBLE'),fliplr(dimid(ids)));	
			netcdf.defVarDeflate(ncid,varId,true,true,Compression);
			if length(ids)>1
				netcdf.putVar(ncid,varId,permute(Obj.(var),length(ids):-1:1)); % we need to reverse the dimension order because Matlab netcdf API saves data in the reverse order
			else
				netcdf.putVar(ncid,varId,Obj.(var));
			end
			for jj=1:length(f)
				if ~isempty(strfind(f{jj},[var '_']))
					netcdf.putAtt(ncid,varId,f{jj}(strfind(f{jj},[var '_'])+length([var '_']):end),Obj.(f{jj}));
				end
			end		
		end
	end

%% Save data variables and their attributes
fd=fieldnames(Obj.Dimensions.Data);
fod=fieldnames(Obj.Data);
	
	for ii=1:length(fd)
		var=fd{ii};
		if isempty(strfind(var,'_'))	% skip all attributes				
			ids=cell2mat(regexp(Dims,cellstr((Obj.Dimensions.Data.(var))')));
			varId = netcdf.defVar(ncid,['Data.' var],netcdf.getConstant('NC_DOUBLE'),fliplr(dimid(ids)));	
			netcdf.defVarDeflate(ncid,varId,true,true,Compression);
			if length(ids)>1
				netcdf.putVar(ncid,varId,permute(Obj.Data.(var),length(ids):-1:1)); % we need to reverse the dimension order because Matlab netcdf API saves data in the reverse order
			else
				netcdf.putVar(ncid,varId,Obj.Data.(var));
			end
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
	for ii=1:length(ME.stack)
		disp(ME.stack(ii));
	end
	error(['Error processing ' var 10 ...
					'Error message: ' ME.message]);

end
netcdf.close(ncid);
	
