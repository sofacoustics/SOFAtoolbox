function NETCDFsave(filename,Obj,Compression)
%NETCDFSAVE
%   NETCDFsave(filename,Dataset,Compression) saves all data and metadata to
%   a SOFA file.

% SOFA API - function netcdf/NETCDFsave
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

% Piotr Majdak, 9.4.2013

% If we are running octave we have to import the NETCDF namespace, in order to
% run functions like netcdf.getConstant
if exist('OCTAVE_VERSION','builtin')
    import_netcdf;
end

%% Global definitions
glob='GLOBAL_';
% Dims='IRENMCQ'; % dimensions

globid=netcdf.getConstant('GLOBAL');

try 
	var='file creation';
  mode = netcdf.getConstant('NETCDF4');
%   mode = netcdf.getConstant('clobber');
%   mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
	ncid = netcdf.create(filename,mode);

%% Save global attributes
  f=fieldnames(Obj);

	for ii=1:length(f)
		if ~isempty(strfind(f{ii},glob))
			var=f{ii};
      netcdf.putAtt(ncid,globid,var(strfind(var,glob)+length(glob):end),Obj.(var));
		end
	end
	
%% Define dimensions

  Dims=cell2mat(fieldnames(rmfield(Obj.API,'Dimensions'))');
	dimid=nan(size(Dims));
	dimsize=nan(size(Dims));
	for ii=1:length(Dims)
		var=Dims(ii);
		dimid(ii) = netcdf.defDim(ncid,Dims(ii),Obj.API.(var)); 
		dimsize(ii)=Obj.API.(var);
  end
  Sdim=strfind(Dims,'S'); % find the index with string dimension
  if isempty(Sdim), Sdim=-1, end; % mark as -1 if not existing

%% Define metadata variables and their attributes
Dimensions=rmfield(Obj.API.Dimensions,'Data');
fv=fieldnames(Dimensions);

	for ii=1:length(fv)
		var=fv{ii};
		if isempty(strfind(var,'_')) % skip all attributes
			ids=cell2mat(regexp(Dims,cellstr((Dimensions.(var))')));
      if find(ids==Sdim) % array of strings or numerics?
        varId(ii) = netcdf.defVar(ncid,var,netcdf.getConstant('NC_CHAR'),fliplr(dimid(ids)));	
      else
        varId(ii) = netcdf.defVar(ncid,var,netcdf.getConstant('NC_DOUBLE'),fliplr(dimid(ids)));	
      end
			netcdf.defVarDeflate(ncid,varId(ii),true,true,Compression);
			for jj=1:length(f)
				if ~isempty(strfind(f{jj},[var '_']))
					netcdf.putAtt(ncid,varId(ii),f{jj}(strfind(f{jj},[var '_'])+length([var '_']):end),Obj.(f{jj}));
				end
			end		
		end
	end

%% Define data variables and their attributes
fd=fieldnames(Obj.API.Dimensions.Data);
fod=fieldnames(Obj.Data);
	
	for ii=1:length(fd)
		var=fd{ii};
		if isempty(strfind(var,'_'))	% skip all attributes				
			ids=cell2mat(regexp(Dims,cellstr((Obj.API.Dimensions.Data.(var))')));
      if find(ids==Sdim) % array of strings or numerics?
        varIdD(ii) = netcdf.defVar(ncid,['Data.' var],netcdf.getConstant('NC_CHAR'),fliplr(dimid(ids)));	
      else
        varIdD(ii) = netcdf.defVar(ncid,['Data.' var],netcdf.getConstant('NC_DOUBLE'),fliplr(dimid(ids)));	
      end
			netcdf.defVarDeflate(ncid,varIdD(ii),true,true,Compression);
			for jj=1:length(fod)
				if ~isempty(strfind(fod{jj},[var '_']))
					netcdf.putAtt(ncid,varIdD(ii),fod{jj}(strfind(fod{jj},[var '_'])+length([var '_']):end),Obj.Data.(fod{jj}));
				end
			end		
		end
  end

%% End of definition
	netcdf.endDef(ncid);
  
%% Save metadata variables
Dimensions=rmfield(Obj.API.Dimensions,'Data');
fv=fieldnames(Dimensions);

	for ii=1:length(fv)
		var=fv{ii};
		if isempty(strfind(var,'_')) % skip all attributes
			ids=cell2mat(regexp(Dims,cellstr((Dimensions.(var))')));
			if length(ids)>1
        if iscell(Obj.(var))
          c=char(permute(Obj.(var),length(ids):-1:1));
          s=size(c);
          ext=zeros([s(1:end-1) Obj.API.S-s(end)]);
          netcdf.putVar(ncid,varId(ii),[c ext]);
        else  
          netcdf.putVar(ncid,varId(ii),permute(Obj.(var),length(ids):-1:1)); % we need to reverse the dimension order because Matlab netcdf API saves data in the reverse order
        end
      else
        if iscell(Obj.(var))
          c=char(Obj.(var));  % convert to character array
          s=size(c);      
          ext=zeros([s(1:end-1) Obj.API.S-s(end)]); % extend along the S dimension
          netcdf.putVar(ncid,varId(ii),[c ext]);  
        else
          netcdf.putVar(ncid,varId(ii),Obj.(var));
        end
			end
		end
	end

%% Save data variables
fd=fieldnames(Obj.API.Dimensions.Data);
fod=fieldnames(Obj.Data);
	
	for ii=1:length(fd)
		var=fd{ii};
		if isempty(strfind(var,'_'))	% skip all attributes				
			ids=cell2mat(regexp(Dims,cellstr((Obj.API.Dimensions.Data.(var))')));
			if length(ids)>1
        if iscell(Obj.Data.(var))
          c=char(permute(Obj.Data.(var),length(ids):-1:1)); % we need to reverse the dimension order because Matlab netcdf API saves data in the reverse order
          s=size(c);
          ext=zeros([s(1:end-1) Obj.API.S-s(end)]);
          netcdf.putVar(ncid,varIdD(ii),[c ext]); 
        else
          netcdf.putVar(ncid,varIdD(ii),permute(Obj.Data.(var),length(ids):-1:1)); % we need to reverse the dimension order because Matlab netcdf API saves data in the reverse order
        end
      else
        if iscell(Obj.Data.(var))
          c=char(Obj.Data.(var));  % convert to character array
          s=size(c);      
          ext=zeros([s(1:end-1) Obj.API.S-s(end)]); % extend along the S dimension          
          netcdf.putVar(ncid,varIdD(ii),[c ext]);
        else
          netcdf.putVar(ncid,varIdD(ii),Obj.Data.(var));
        end
			end
		end
	end
  
catch ME
% 	if ~strcmp(ME.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
		netcdf.close(ncid);
% 	end
	for ii=1:length(ME.stack)
		disp(ME.stack(ii));
	end
	error(['Error processing ' var 10 ...
					'Error message: ' ME.message]);

end
netcdf.close(ncid);
	
