function Obj = SOFAupdateDimensions(Obj)
%SOFAupdateDimensions
%   Obj = SOFAupdateDimensions(Obj) updates the dimensions in the SOFA
%   structure
%
%   Obj is a struct containing the data and meta.
%		The dimension sizes are created as .API.DimSize and updated corresponding to the
%		conventions

% SOFA API - function SOFAupdateDimensions
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Get conventions with allowed dimensions
OC = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'a');

%% Update dimension sizes

  % fix dimension sizes
Obj.API.DimSize.I=1;
Obj.API.DimSize.C=3;
  % variable-dependent dimension sizes
dims='renm'; 
  % check all metadata variables
f=fieldnames(rmfield(OC.API.Dimensions,'Data'));
for ii=1:length(dims)
	for jj=1:length(f)
		dim=strfind(OC.API.Dimensions.(f{jj}),dims(ii));
		if iscell(dim), dim=cell2mat(dim); end;
		if ~isempty(dim)
			Obj.API.DimSize.(upper(dims(ii)))=size(Obj.(f{jj}),dim(1));
			break;
		end
	end
end
  % check all data variables
fd=fieldnames(OC.API.Dimensions.Data);
for ii=1:length(dims)
	for jj=1:length(fd)
		dim=strfind(OC.API.Dimensions.Data.(fd{jj}),dims(ii));
		if iscell(dim), dim=cell2mat(dim); end;
		if ~isempty(dim)
			Obj.API.DimSize.(upper(dims(ii)))=size(Obj.Data.(fd{jj}),dim(1));
			break;
		end
	end
end

%% Update the dimensions of metadata variables
X=rmfield(Obj,{'Data','API'});
if isfield(X,'PRIVATE'), X=rmfield(X,'PRIVATE'); end
Xf=fieldnames(X);
for ii=1:length(Xf)
	if isempty(strfind(Xf{ii},'_')),	% is not an attribute...
		if isfield(OC.API.Dimensions, Xf{ii}), % is a known variable		
%       disp(Xf{ii});
			dim=OC.API.Dimensions.(Xf{ii});
			if ~iscell(dim), dim={dim}; end;
			dim=checkdim(Obj,dim,size(Obj.(Xf{ii})));
			if isempty(dim),
				error([Xf{ii} ': dimension could not be matched.']);
			else
				Obj.API.Dimensions.(Xf{ii})=dim;
			end
		else % is a user-defined variable						
			if ~isfield(Obj.API.Dimensions,Xf{ii}),
				error([Xf{ii} ' seems to be a user-defined variable without a dimension.']);
      else
        dim=Obj.API.Dimensions.(Xf{ii});
        dim=checkdim(Obj,{dim},size(Obj.(Xf{ii})));
        if isempty(dim),
          error([Xf{ii} ': dimension does not match.']);
        end
			end
		end		
	end
end

%% Update the dimensions of data variables
Xf=fieldnames(Obj.Data);
for ii=1:length(Xf)
	if isempty(strfind(Xf{ii},'_')),	% is not an attribute...
		if isfield(OC.API.Dimensions.Data, Xf{ii}), 			% is a known variable
			dim=OC.API.Dimensions.Data.(Xf{ii}); 
			if ~iscell(dim), dim={dim}; end;
			dim=checkdim(Obj,dim,size(Obj.Data.(Xf{ii})));
			if isempty(dim),
				error(['Data.' Xf{ii} ': dimension could not be matched.']);
			else
				Obj.API.Dimensions.Data.(Xf{ii})=dim;
			end
		else
			error(['Unknown data variable' Xf{ii} '.']);
		end		
	end
end

%% Get the sizes of the dimension variables according the dimension variables in str
function vec=getdim(Obj,str)
	vec=arrayfun(@(f)(Obj.(f)),upper(str));

%% dims is a cell array with allowed dimensions. dimA is a vector with the actual dimensions.
% dim is a string with the matching dimension
function dim=checkdim(Obj,dims,dimA)
dim=[];
for jj=1:length(dims)
  dimS=dims{jj};
  if length(dimS)==1, dimS=[dimS 'I']; end; % 1D required, but Matlab is always 2D at least.
	dimR=getdim(Obj.API.DimSize,dimS);
	if length(dimA)==length(dimR), % the same size?
		if dimA==dimR, dim=upper(dims{jj}); break; end;	% found!
	elseif length(dimA)<length(dimR)	% extend the size?
		if [dimA ones(1,length(dimR)-length(dimA))]==dimR, dim=upper(dims{jj}); break; end; % found!
	end
end