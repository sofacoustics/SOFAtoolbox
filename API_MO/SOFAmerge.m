function [C, log] = SOFAmerge(A,B)
%SOFAmerge
%   [C, log] = SOFAmerge(A,B) merges the SOFA objects A and B to a single one, C.
%
%   A and B are structs containing the data and meta. A and B
%   must be of the same SOFA conventions.
% 
%   C is a struct containing merged data
%   log contains log data

% SOFA API - function SOFAupdateDimensions
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Initial check 
if ~strcmp(A.GLOBAL_SOFAConventions,B.GLOBAL_SOFAConventions)
	error('Both SOFA objects must use the same SOFA conventions');
end
if A.API.N~=B.API.N
	error('Data size of both SOFA objects must be the same');
end
Def = SOFAdefinitions;
log={''};
OC = SOFAgetConventions(A.GLOBAL_SOFAConventions,'a');

%% Update dimensions
A=SOFAupdateDimensions(A);
B=SOFAupdateDimensions(B);

%% create field names which have to be checked
C=rmfield(B,{'API','Data'});
if isfield(C,'PRIVATE'), C=rmfield(C,'PRIVATE'); end
Bf=fieldnames(C);

%% Copy and merge metadata variables
C=A;
for ii=1:size(Bf,1)
	if ~isfield(A,Bf{ii})
		C.(Bf{ii}) = B.(Bf{ii});	% field in B but not in A. Simple copy.		
	else	% here we have a potential conflict and have to merge
		if strfind(Bf{ii},'_')	% is it an attribute?
			if strcmp(A.(Bf{ii}),B.(Bf{ii})), 
				C.(Bf{ii}) = B.(Bf{ii});	% content the same, no conflict
			else
				switch Bf{ii}
					case 'GLOBAL_TimeCreated'	% use the oldest date
						dateNew=datenum(A.GLOBAL_TimeCreated,Def.dateFormat);
						if datenum(B.GLOBAL_TimeCreated,Def.dateFormat)<dateNew, dateNew=datenum(B.GLOBAL_TimeCreated,Def.dateFormat); end;
						C.(Bf{ii}) = datestr(dateNew,Def.dateFormat);
						log{end+1}=[Bf{ii} ' set to ' C.(Bf{ii})];
					case 'GLOBAL_TimeModified' % now
						C.(Bf{ii}) = datestr(now,Def.dateFormat);
						log{end+1}=[Bf{ii} ' updated'];
					otherwise
						C.(Bf{ii}) = [A.(Bf{ii}) '; ' B.(Bf{ii})]; % concatate [A; B]
						log{end+1}=[Bf{ii} ' merged'];
				end
			end
		else	% a variable
			if isfield(OC.API.Dimensions, Bf{ii})	% is a known variable?
				AExp=SOFAexpand(A,Bf{ii});
				BExp=SOFAexpand(B,Bf{ii});
				dim=strfind(AExp.API.Dimensions.(Bf{ii}),'M');	
				if isempty(dim),
					error([Bf{ii} ' can not be merged because it does not depend on M']);
				end
				C.(Bf{ii})=cat(dim,AExp.(Bf{ii}),BExp.(Bf{ii}));
				log{end+1}=[Bf{ii} ' expanded and merged'];
			else	% user-defined variable, dimensions must be stated
				if ~isfield(A.API.Dimensions, Bf{ii})
					error(['Dimension missing for ' Bf{ii} ' in A.']); end
				if ~isfield(B.API.Dimensions, Bf{ii})
					error(['Dimension missing for ' Bf{ii} ' in B.']); end
				dim=strfind(A.API.Dimensions.(Bf{ii}),'M');	
				if ~isempty(dim),
          C.(Bf{ii})=cat(dim,A.(Bf{ii}),B.(Bf{ii})); % depends on M, merge
          log{end+1}=[Bf{ii} ' merged'];
        else  % not M-dependend, must be identical in A and B
          if prod(A.(Bf{ii})==B.(Bf{ii}))==1,
            C.(Bf{ii})=A.(Bf{ii});
            log{end+1}=[Bf{ii} ' identical'];
          else
            error([Bf{ii} ' must depend on M or be equal in both objects.']);
          end
        end
			end
		end
	end
end

%% Copy and merge Data variables
Bf=fieldnames(B.Data);
for ii=1:size(Bf,1)
	if ~isfield(A.Data,Bf{ii})
		C.Data.(Bf{ii}) = B.Data.(Bf{ii});	% field in B but not in A. Simple copy.		
	else	% here we have a potential conflict and have to merge
		if strfind(Bf{ii},'_')	% is it an attribute?
			if strcmp(A.Data.(Bf{ii}),B.Data.(Bf{ii})), 
				C.Data.(Bf{ii}) = B.Data.(Bf{ii});	% content the same, no conflict
			else
				C.Data.(Bf{ii}) = [A.Data.(Bf{ii}) '; ' B.Data.(Bf{ii})]; % concatate [A; B]
				log{end+1}=['Data.' Bf{ii} ' merged'];
			end
		else	% a variable in Data
			if isfield(OC.API.Dimensions.Data,Bf{ii})	% is a known variable?
				dim=strfind(A.API.Dimensions.Data.(Bf{ii}),'M'); % is a matrix
				if ~isempty(dim),
          C.Data.(Bf{ii})=cat(dim,A.Data.(Bf{ii}),B.Data.(Bf{ii})); % depends on M, cat
          log{end+1}=['Data.' Bf{ii} ' merged'];
        else  % not M-dependend, must be identical in A and B
          if all(A.Data.(Bf{ii})==B.Data.(Bf{ii}))==1,
            C.Data.(Bf{ii})=A.Data.(Bf{ii});
            log{end+1}=['Data.' Bf{ii} ' identical'];
          else
            error(['Data.' Bf{ii} ' must depend on M or be equal in both objects.']);
          end
        end
			else	% user-defined variable, dimensions must be stated
				if ~isfield(A.API.Dimensions.Data, Bf{ii})
					error(['Dimension missing for Data.' Bf{ii} ' in A.']); end
				if ~isfield(B.API.Dimensions.Data, Bf{ii})
					error(['Dimension missing for Data.' Bf{ii} ' in B.']); end
				dim=strfind(A.API.Dimensions.Data.(Bf{ii}),'M');	
				if ~isempty(dim),
          C.Data.(Bf{ii})=cat(dim,A.Data.(Bf{ii}),B.Data.(Bf{ii})); % depends on M, cat
          log{end+1}=['Data.' Bf{ii} ' merged'];
        else  % not M-dependend, must be identical in A and B
          if prod(A.Data.(Bf{ii})==B.Data.(Bf{ii}))==1,
            C.Data.(Bf{ii})=A.Data.(Bf{ii});
            log{end+1}=['Data.' Bf{ii} ' identical'];
          else
            error(['Data.' Bf{ii} ' must depend on M or be equal in both objects.']);
          end
        end
			end
		end
	end
end

%% Update the new dimensions and finish
C=SOFAcompact(C);
C=SOFAupdateDimensions(C);
if length(log)>1, log=log(2:end); else log={}; end;

%% Get the sizes of the dimension variables according the dimension variables in str
function vec=getdim(Obj,str)
	vec=arrayfun(@(f)(Obj.API.(f)),upper(str));
