function [C, log] = SOFAmerge(A,B)
%SOFAmerge
%   C = SOFAmerge(A,B) merges the SOFA objects A and B to a single one, C.
%
%   A and B are structs containing the data and meta. A and B
%   must be of the same SOFA conventions.

% SOFA API - function SOFAupdateDimensions
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% Initial check 
if ~strcmp(A.GLOBAL_SOFAConventions,B.GLOBAL_SOFAConventions)
	error('Both SOFA objects must use the same SOFA conventions');
end
if A.N~=B.N
	error('Data size of both SOFA objects must be the same');
end
log={''};
OC = SOFAgetConventions(A.GLOBAL_SOFAConventions,'a');

%% Update dimensions
A=SOFAupdateDimensions(A);
B=SOFAupdateDimensions(B);

%% create field names which have to be checked
If={'Dimensions','Data','I','R','E','N','M','C','Q'};
C=B;
for ii=1:length(If)
	if isfield(C,If{ii}), C=rmfield(C,If{ii}); end
end
Bf=fieldnames(C);

%% Copy and merge w/o data
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
					case 'GLOBAL_DatabaseTimeCreated'	% use the oldest date
						dateNew=datenum(A.GLOBAL_DatabaseTimeCreated);
						if datenum(B.GLOBAL_DatabaseTimeCreated)<dateNew, dateNew=datenum(B.GLOBAL_DatabaseTimeCreated); end;
						C.(Bf{ii}) = datestr(dateNew,31);
						log{end+1}=[Bf{ii} ' set to ' C.(Bf{ii})];
					case 'GLOBAL_DatabaseTimeModified' % now
						C.(Bf{ii}) = datestr(now,31);
						log{end+1}=[Bf{ii} ' updated'];
					otherwise
						C.(Bf{ii}) = [A.(Bf{ii}) '; ' B.(Bf{ii})]; % concatate [A; B]
						log{end+1}=[Bf{ii} ' merged'];
				end
			end
		else	% a variable
			if isfield(OC.Dimensions, Bf{ii})	% is a known variable?
				AExp=SOFAexpand(A,Bf{ii});
				BExp=SOFAexpand(B,Bf{ii});
				dim=strfind(AExp.Dimensions.(Bf{ii}),'M');	
				if isempty(dim),
					error([Bf{ii} ' can not be merged because it does not depend on M']);
				end
				C.(Bf{ii})=cat(dim,AExp.(Bf{ii}),BExp.(Bf{ii}));
				log{end+1}=[Bf{ii} ' expanded and merged'];
			else	% user-defined variable, dimensions must be stated
				if ~isfield(A.Dimensions, Bf{ii})
					error(['Dimension missing for ' Bf{ii} ' in A.']); end
				if ~isfield(B.Dimensions, Bf{ii})
					error(['Dimension missing for ' Bf{ii} ' in B.']); end
				dim=strfind(A.Dimensions.(Bf{ii}),'M');	
				if isempty(dim),
					error([Bf{ii} ' can not be merged because it does not depend on M']);
				end
				C.(Bf{ii})=cat(dim,A.(Bf{ii}),B.(Bf{ii}));
				log{end+1}=[Bf{ii} ' expanded and merged'];
			end
		end
	end
end

%% Copy and merge data
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
			if isfield(OC.Dimensions.Data,Bf{ii})	% is a known variable?
				if strcmp(A.Dimensions.Data.(Bf{ii}),'I'),	% must be scalar?
					if A.Data.(Bf{ii})~=B.Data.(Bf{ii}),
						error(['Data.' Bf{ii} ' must be scalar and is not the same in A and B']); end
					C.Data.(Bf{ii})=A.Data.(Bf{ii});
					continue;
				end
				dim=strfind(A.Dimensions.Data.(Bf{ii}),'M'); % is a matrix
				if isempty(dim),
					error(['Data.' Bf{ii} ' can not be merged because it does not depend on M']); end
				C.Data.(Bf{ii})=cat(dim,A.Data.(Bf{ii}),B.Data.(Bf{ii}));
				log{end+1}=['Data.' Bf{ii} ' merged'];
			else	% user-defined variable, dimensions must be stated
				if ~isfield(A.Dimensions.Data, Bf{ii})
					error(['Dimension missing for Data.' Bf{ii} ' in A.']); end
				if ~isfield(B.Dimensions.Data, Bf{ii})
					error(['Dimension missing for Data.' Bf{ii} ' in B.']); end
				dim=strfind(A.Dimensions.Data.(Bf{ii}),'M');	
				if isempty(dim),
					error(['Data.' Bf{ii} ' can not be merged because it does not depend on M']); end
				C.Data.(Bf{ii})=cat(dim,A.Data.(Bf{ii}),B.Data.(Bf{ii}));
				log{end+1}=['Data.' Bf{ii} ' merged'];
			end
		end
	end
end

%% Update the new dimensions and finish
C=SOFAupdateDimensions(C);
if length(log)>1, log=log(2:end); else log={}; end;

%% Get the sizes of the dimension variables according the dimension variables in str
function vec=getdim(Obj,str)
	vec=arrayfun(@(f)(Obj.(f)),upper(str));