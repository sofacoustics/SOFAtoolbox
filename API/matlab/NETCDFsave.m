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

% Last Update: Michael Mihocic, 04.2013
% Piotr Majdak, 9.4.2013



%% Global definitions
glob='GLOBAL_';
dims='MRNECQ';
f=fieldnames(Obj);
fv=fieldnames(Var);
fd=fieldnames(DataVar);

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
% 	f=fieldnames(DataVar);
% 	for ii=1:length(f)
% 		if isempty(strfind(f{ii},'_'))	% skip all attributes
% 			SaveDataVariableAndAttributes(ncid,Obj.Data,f{ii},Cid(DataVar.(f{ii})));
% 		end
% 	end
	
catch ME
	if ~strcmp(ME.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
		netcdf.close(ncid);
	end
	error([ME.message ' Error in line ' num2str(ME.stack.line)]);
end
	netcdf.close(ncid);
	
	
function SaveVariableAndAttributes(ncid,Obj,var,dimid)
disp(['VAR:' var num2str(size(Obj.(var)))]);
VarId = netcdf.defVar(ncid,var,netcdf.getConstant('NC_FLOAT'),dimid);	
netcdf.putVar(ncid,VarId,Obj.(var));
f=fieldnames(Obj);
for ii=1:length(f)
	if ~isempty(strfind(f{ii},[var '_']))
		netcdf.putAtt(ncid,VarId,f{ii}(strfind(f{ii},[var '_'])+length([var '_']):end),Obj.(f{ii}));
	end
end

function SaveDataVariableAndAttributes(ncid,Obj,var,dimid)
disp(['DATA: ' 'Data.' var]);
VarId = netcdf.defVar(ncid,['Data.' var],netcdf.getConstant('NC_FLOAT'),dimid);	
netcdf.putVar(ncid,VarId,Obj.(var));
f=fieldnames(Obj);
for ii=1:length(f)
	if ~isempty(strfind(f{ii},[var '_']))
		netcdf.putAtt(ncid,VarId,f{ii}(strfind(f{ii},[var '_'])+length([var '_']):end),Obj.(f{ii}));
	end
end



%         % define some constants and fixed dimensions
%         float = netcdf.getConstant('NC_FLOAT');
%         ScalarDimId = netcdf.defDim(ncid,dimNames.Scalar,1);
%         UnlimitedDimId= netcdf.defDim(ncid,dimNames.Unlimited,netcdf.getConstant('NC_UNLIMITED'));

%         
%         
%         for ii=1:numVars % loop through all input variables
%             currentVarName = varNames{ii};
% %             disp(currentVarName);
%             currentVarValue = Obj.(varNames{ii});
%             
%             switch currentVarName
%               
% 
%                 case 'RLongName'
%                     VarId = netcdf.inqVarID(ncid,'R');
%                     netcdf.putAtt(ncid,VarId,'LongName',currentVarValue);
%                 case 'ELongName'
%                     VarId = netcdf.inqVarID(ncid,'E');
%                     netcdf.putAtt(ncid,VarId,'LongName',currentVarValue);
%                 case 'MLongName'
%                     VarId = netcdf.inqVarID(ncid,'M');
%                     netcdf.putAtt(ncid,VarId,'LongName',currentVarValue);
%                 case 'CLongName'
%                     VarId = netcdf.inqVarID(ncid,'C');
%                     netcdf.putAtt(ncid,VarId,'LongName',currentVarValue);                    
%                 case 'NLongName'
%                     VarId = netcdf.inqVarID(ncid,'N');
%                     netcdf.putAtt(ncid,VarId,'LongName',currentVarValue);                
%                 case 'NUnits'
%                     VarId = netcdf.inqVarID(ncid,'N');
%                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue);  
% %                 case 'SamplingRateUnits'
% %                     VarId = netcdf.inqVarID(ncid,'Data.SamplingRate');
% %                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue);                    
%                 case 'ListenerPositionType'
%                     VarId = netcdf.inqVarID(ncid,'ListenerPosition');
%                     netcdf.putAtt(ncid,VarId,'Type',currentVarValue);   
%                 case 'ListenerPositionUnits'
%                     VarId = netcdf.inqVarID(ncid,'ListenerPosition');
%                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue);  
%                 case 'ListenerUpType'
%                     VarId = netcdf.inqVarID(ncid,'ListenerUp');
%                     netcdf.putAtt(ncid,VarId,'Type',currentVarValue);   
%                 case 'ListenerUpUnits'
%                     VarId = netcdf.inqVarID(ncid,'ListenerUp');
%                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue);   
%                 case 'ListenerViewType'
%                     VarId = netcdf.inqVarID(ncid,'ListenerView');
%                     netcdf.putAtt(ncid,VarId,'Type',currentVarValue);   
%                 case 'ListenerViewUnits'
%                     VarId = netcdf.inqVarID(ncid,'ListenerView');
%                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue);   
%                 case 'ListenerRotationType'
%                     VarId = netcdf.inqVarID(ncid,'ListenerRotation');
%                     netcdf.putAtt(ncid,VarId,'Type',currentVarValue);   
%                 case 'ListenerRotationUnits'
%                     VarId = netcdf.inqVarID(ncid,'ListenerRotation');
%                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue);      
%                 case 'ReceiverPositionType'
%                     VarId = netcdf.inqVarID(ncid,'ReceiverPosition');
%                     netcdf.putAtt(ncid,VarId,'Type',currentVarValue);   
%                 case 'ReceiverPositionUnits'
%                     VarId = netcdf.inqVarID(ncid,'ReceiverPosition');
%                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue); 
%                 case 'SourcePositionType'
%                     VarId = netcdf.inqVarID(ncid,'SourcePosition');
%                     netcdf.putAtt(ncid,VarId,'Type',currentVarValue);   
%                 case 'SourcePositionUnits'
%                     VarId = netcdf.inqVarID(ncid,'SourcePosition');
%                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue);      
%                 case 'SourceUpType'
%                     VarId = netcdf.inqVarID(ncid,'SourceUp');
%                     netcdf.putAtt(ncid,VarId,'Type',currentVarValue);   
%                 case 'SourceUpUnits'
%                     VarId = netcdf.inqVarID(ncid,'SourceUp');
%                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue);          
%                 case 'SourceViewType'
%                     VarId = netcdf.inqVarID(ncid,'SourceView');
%                     netcdf.putAtt(ncid,VarId,'Type',currentVarValue);   
%                 case 'SourceViewUnits'
%                     VarId = netcdf.inqVarID(ncid,'SourceView');
%                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue);      
%                 case 'EmitterPositionType'
%                     VarId = netcdf.inqVarID(ncid,'EmitterPosition');
%                     netcdf.putAtt(ncid,VarId,'Type',currentVarValue);   
%                 case 'EmitterPositionUnits'
%                     VarId = netcdf.inqVarID(ncid,'EmitterPosition');
%                     netcdf.putAtt(ncid,VarId,'Units',currentVarValue);                     
%                 case 'RoomDAEFileNameDescription'
%                     VarId = netcdf.inqVarID(ncid,'RoomDAEFileName');
%                     netcdf.putAtt(ncid,VarId,'Description',currentVarValue); 
%                     
%                 otherwise
% %                     variables                    
%             % ------------------ check and prepare variables ------------------
%             % -- convert all strings to cells --
%             if ~isnumeric(currentVarValue) && ~isstruct(currentVarValue) % if currentVarValue is a string
%                 currentVarValue = cellstr(currentVarValue);
%             end
%             % dimensions (for length of string) if a cell only contains one string
%             if ~isnumeric(currentVarValue) && ~isstruct(currentVarValue) % -- if currentVarValue is a string
%                 if size(currentVarValue,1) == 1 && size(currentVarValue,2) == 1 % [1 1]
%                     DimId = netcdf.defDim(ncid,[currentVarName 'DIM'],length(currentVarValue{1}));
%                 end
%                 if size(currentVarValue,1) == 1 && size(currentVarValue,2) > 1 % [1 x]
%                     xDimId = netcdf.defDim(ncid,[currentVarName 'xDIM'],size(currentVarValue,2));
%                     lengths=zeros(size(currentVarValue,2));
%                     for n=1:size(currentVarValue,2) % go through all strings up to x
%                         lengths(n) = length(currentVarValue{n}); % store all string lengths
%                     end
%                     % length of dimension is maximum of all string lengths
%                     DimId = netcdf.defDim(ncid,[currentVarName 'DIM'],max(lengths));
%                 end
%                 if size(currentVarValue,1) == M && size(currentVarValue,2) == 1 % [M 1]
%                     lengths=zeros(M);
%                     for n=1:M % go through all strings up to x
%                         lengths(n) = length(currentVarValue{n});
%                     end
%                     DimId = netcdf.defDim(ncid,[currentVarName 'DIM'],max(lengths));
%                 end
%                 if size(currentVarValue,1) == M && size(currentVarValue,2) > 1 % [M x]
%                     xDimId = netcdf.defDim(ncid,[currentVarName 'xDIM'],size(currentVarValue,2));
%                     lengths=zeros(M,size(currentVarValue,2));
%                     for n=1:M % go through all strings up to M and x (2D)
%                         for m=1:size(currentVarValue,2)
%                             lengths(n,m) = length(currentVarValue{n,m});
%                         end
%                     end
%                     DimId = netcdf.defDim(ncid,[currentVarName 'DIM'],max(max(lengths)));
%                 end
%             % dimensions of length x for normal, numeric variables, [1 x] or [M x]
%             elseif ~(strcmp(currentVarName,'Data') || sum(strcmp(currentVarName,SourceListenerVars)) || ...
%                      sum(strcmp(currentVarName,TransmitterReceiverVars)))
%                 if size(currentVarValue,2) > 1
%                     DimId = netcdf.defDim(ncid,[currentVarName 'xDIM'],size(currentVarValue,2)); 
%                 end
%             end
% 
%             % ------------------------ variables ---------------------------
%             if ~isnumeric(currentVarValue) && ~isstruct(currentVarValue) % --- define string variables ---
%                 % string variable, [1 1], [1 x], [M 1], [M x]
%                 if size(currentVarValue,1) == 1 && size(currentVarValue,2) == 1 % [1 1]
%                     VarId = netcdf.defVar(ncid,currentVarName,2,[ScalarDimId DimId]);
%                 end
%                 if size(currentVarValue,1) == 1 && size(currentVarValue,2) > 1 % [1 x]
%                     VarId = netcdf.defVar(ncid,currentVarName,2,[ScalarDimId xDimId DimId]);
%                 end
%                 if size(currentVarValue,1) == M && size(currentVarValue,2) == 1 % [M 1]
%                     VarId = netcdf.defVar(ncid,currentVarName,2,[MDimId ScalarDimId DimId]);
%                 end
%                 if size(currentVarValue,1) == M && size(currentVarValue,2) > 1 % [M x]
%                     VarId = netcdf.defVar(ncid,currentVarName,2,[MDimId xDimId DimId]);
%                 end
%             else % --- define numeric variables ---
%                 if strcmp(currentVarName,'Data') % -- Data, float, [M R N]
%                     dataTypes=fieldnames(currentVarValue);
%                     for jj=1:length(dataTypes)
%                         switch dataTypes{jj}
%                             case 'SamplingRateUnits'
%                                 VarId = netcdf.inqVarID(ncid,'Data.SamplingRate');
%                                 netcdf.putAtt(ncid,VarId,'Units',currentVarValue.SamplingRateUnits);   
%                             case 'SamplingRate'
%                                 VarId=netcdf.defVar(ncid,['Data.' dataTypes{jj}],float,ScalarDimId);
% %                                 VarId(jj)=netcdf.defVar(ncid,['Data.' dataTypes{jj}],'double',[MDimId RDimId NDimId]);
%                                 netcdf.defVarDeflate(ncid,VarId,true,true,Compression); % Compression of data
%                                 netcdf.putVar(ncid,VarId,currentVarValue.(dataTypes{jj}));
%                                 issaved=1;                                
%                             otherwise
%                                 VarId=netcdf.defVar(ncid,['Data.' dataTypes{jj}],'double',[MDimId RDimId NDimId]);
% %                                 VarId(jj)=netcdf.defVar(ncid,['Data.' dataTypes{jj}],'double',[MDimId RDimId NDimId]);
%                                 netcdf.defVarDeflate(ncid,VarId,true,true,Compression); % Compression of data
%                                 netcdf.putVar(ncid,VarId,currentVarValue.(dataTypes{jj}));
%                                 issaved=1;
%                         end
%                     end 
%                 elseif sum(strcmp(currentVarName,SourceListenerVars))
%                     % -- positions and vectors, float, [1 3] or [M 3]
%                     if size(currentVarValue,1) > 1
%                         VarId = netcdf.defVar(ncid,currentVarName,float,[MDimId CoordDimId]);
%                     else
%                         VarId = netcdf.defVar(ncid,currentVarName,float,[ScalarDimId CoordDimId]);
%                     end
% 
%                 elseif sum(strcmp(currentVarName,'ReceiverPosition'))
%                     % receiver position, float, [R 3], [R 3 M]
%                     if (size(currentVarValue,3) > 1)
%                         VarId = netcdf.defVar(ncid,currentVarName,float,[RDimId CoordDimId MDimId]);
%                                     else
%                         VarId = netcdf.defVar(ncid,currentVarName,float,[RDimId CoordDimId ScalarDimId]);
%                     end
%                 elseif sum(strcmp(currentVarName,'TransmitterPosition'))
%                     % transmitter position, float, [T 3], [T 3 M]
%                     if (size(currentVarValue,3) > 1)
%                         VarId = netcdf.defVar(ncid,currentVarName,float,[TDimId CoordDimId MDimId]);
%                                     else
%                         VarId = netcdf.defVar(ncid,currentVarName,float,[TDimId CoordDimId ScalarDimId]);
%                     end
%                 else % "normal" numeric variables, float, [1 1], [M 1], [1 x], [M x]
%                     if (size(currentVarValue,1) == 1) && (size(currentVarValue,2) == 1)
%                         VarId = netcdf.defVar(ncid,currentVarName,float,ScalarDimId);
%                     end
%                     if (size(currentVarValue,1) == M) && (size(currentVarValue,2) == 1)
%                         VarId = netcdf.defVar(ncid,currentVarName,float,[MDimId ScalarDimId]);
%                     end
%                     if (size(currentVarValue,1) == 1) && (size(currentVarValue,2) > 1)
%                         VarId = netcdf.defVar(ncid,currentVarName,float,[ScalarDimId DimId]);
%                     end
%                     if (size(currentVarValue,1) == M) && (size(currentVarValue,2) > 1)
%                         VarId = netcdf.defVar(ncid,currentVarName,float,[MDimId DimId]);
%                     end
%                 end
% %                 disp(currentVarName);
%                 if issaved == 0
%                     netcdf.defVarDeflate(ncid,VarId,true,true,Compression); % Compression of data
%                 else
%                     issaved=0;
%                 end
%             end
% 
%             % ------------------- write values to variables -----------------
%             if ~isnumeric(currentVarValue) && ~isstruct(currentVarValue) % write string variables
%                 if size(currentVarValue,1) == 1 && size(currentVarValue,2) == 1 % [1 1]
%                     netcdf.putVar(ncid,VarId,char(currentVarValue));
%                 end
%                 if size(currentVarValue,1) == 1 && size(currentVarValue,2) > 1 % [1 x]
%                     for n=1:size(currentVarValue,2) % write elements of cell to variable one-by-one
%                         netcdf.putVar(ncid,VarId,[0 n-1 0],[1 1 length(currentVarValue{n})],currentVarValue{n});
%                     end
%                 end
%                 if size(currentVarValue,1) == M && size(currentVarValue,2) == 1 % [M 1]
%                     for n=1:M % write elements of cell to variable one-by-one
%                         netcdf.putVar(ncid,VarId,[n-1 0 0],[1 1 length(currentVarValue{n})],currentVarValue{n});
%                     end
%                 end
%                 if size(currentVarValue,1) == M && size(currentVarValue,2) > 1 % [M x]
%                     for n=1:M % write elements of cell to variable one-by-one
%                         for m=1:size(currentVarValue,2)
%                             netcdf.putVar(ncid,VarId,[n-1 m-1 0],[1 1 length(currentVarValue{n,m})],currentVarValue{n,m});
%                         end
%                     end
%                 end
%             elseif strcmp(currentVarName,'Data') % write data variables
% %                 for jj=1:length(dataTypes) % DA IS DA HUND DRIN!!!
% %                     disp(dataTypes{jj});
% %                     netcdf.putVar(ncid,VarId(jj),currentVarValue.(dataTypes{jj}));
% %                 end
%             else % numeric variables
% %                 disp(currentVarValue);
%                 netcdf.putVar(ncid,VarId,currentVarValue);
%             end
%             end
%         end
% 
% 
% 
%           
% 
%             
% %     case 'SimpleDRIRMicArray'
% %         
% %         
% %         
% %     otherwise
% %         
% % end
