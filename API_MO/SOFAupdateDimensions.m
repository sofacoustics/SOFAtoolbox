function Obj = SOFAupdateDimensions(Obj,varargin)
%SOFAupdateDimensions
%   Obj = SOFAupdateDimensions(Obj) updates the dimensions in the SOFA
%   structure
%
%   Obj is a struct containing the data and meta.
%		The dimension sizes are created as .API.X and updated corresponding to the
%		conventions
%   flag is 'nodata' or 'all'; default is 'all'

% 9.8.2014: String support added. 
%
% SOFA API - function SOFAupdateDimensions
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

definput.keyvals.Index=[];
definput.flags.type={'data','nodata'};
[flags,kv]=SOFAarghelper({'Index'},definput,varargin);

%% Get conventions with allowed dimensions
OC = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'a');

%% Add dimensions if required
dims=fieldnames(SOFAdefinitions('dimensions'));
for ii=1:size(dims,1)
  if ~isfield(Obj.API,dims{ii}), Obj.API.(dims{ii})=0; end
end

%% Update dimension sizes from the variables having dominant dimensions sizes
  % fix dimension sizes
Obj.API.I=1;
Obj.API.C=3;
  % check all metadata variables for dominant dimension sizes
dims='renm'; 
f=fieldnames(rmfield(OC.API.Dimensions,'Data'));
for ii=1:length(dims)
	for jj=1:length(f)
		dim=strfind(OC.API.Dimensions.(f{jj}),dims(ii));
		if iscell(dim), dim=cell2mat(dim); end;
		if ~isempty(dim)
			Obj.API.(upper(dims(ii)))=size(Obj.(f{jj}),dim(1));
			break;
		end
	end
end
% check all data variables
if flags.do_data
    fd=fieldnames(OC.API.Dimensions.Data);
    for ii=1:length(dims)
        for jj=1:length(fd)
            dim=strfind(OC.API.Dimensions.Data.(fd{jj}),dims(ii));
            if iscell(dim), dim=cell2mat(dim); end;
            if ~isempty(dim)
                Obj.API.(upper(dims(ii)))=size(Obj.Data.(fd{jj}),dim(1));
                break;
            end
        end
    end
end

%% Update the dimensions of metadata variables
Smax=0;
X=rmfield(Obj,{'Data','API'});
if isfield(X,'PRIVATE'), X=rmfield(X,'PRIVATE'); end
Xf=fieldnames(X);
for ii=1:length(Xf)
	if isempty(strfind(Xf{ii},'_')),	% is not an attribute...
		if isfield(OC.API.Dimensions, Xf{ii}), % is a known variable		
%       disp(Xf{ii});
			dim=OC.API.Dimensions.(Xf{ii});
			if ~iscell(dim), dim={dim}; end;
			[dim,S]=checkdim(Obj,dim,sizecell(Obj.(Xf{ii})));
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
        [dim,S]=checkdim(Obj,{dim},sizecell(Obj.(Xf{ii})));
        if isempty(dim),
          error([Xf{ii} ': dimension does not match.']);
        end        
			end
    end
    Smax=max(Smax,S);
	end
end
%% Update the dimensions of data variables
if flags.do_data
    Xf=fieldnames(Obj.Data);
    for ii=1:length(Xf)
        if isempty(strfind(Xf{ii},'_')),	% is not an attribute...
            if isfield(OC.API.Dimensions.Data, Xf{ii}), 			% is a known variable
                dim=OC.API.Dimensions.Data.(Xf{ii}); 
                if ~iscell(dim), dim={dim}; end;
                [dim,S]=checkdim(Obj,dim,sizecell(Obj.Data.(Xf{ii})));
                if isempty(dim),
                    error(['Data.' Xf{ii} ': dimension could not be matched.']);
                else
                    Obj.API.Dimensions.Data.(Xf{ii})=dim;
                end
                Smax=max(Smax,S);
            else
              if ~isfield(Obj.API.Dimensions.Data,Xf{ii}),
                error([Xf{ii} ' seems to be a user-defined variable without a dimension.']);
              else
                dim=Obj.API.Dimensions.Data.(Xf{ii});
                [dim,S]=checkdim(Obj,{dim},sizecell(Obj.Data.(Xf{ii})));
                if isempty(dim),
                  error(['Data.' Xf{ii} ': dimension does not match.']);
                end  
              end
            end		
        end
    end
end
%% Update the size of the longest string
if Smax>0, Obj.API.S=Smax; end
  
%% Return the size of x. If x is a cell, return the size of the strings in x.
function s=sizecell(x,dim)
if iscell(x)
  s=size(char(x));
  if size(x,1)~=s(1) s=[size(x) s(2)]; end % multidim cellarays: s = [celldim1, celldim2, ... , celldimN, stringdim]
else
  s=size(x);
end

%% Get the sizes of the dimension variables according the dimension variables in str
function vec=getdim(Obj,str)
	vec=arrayfun(@(f)(Obj.(f)),upper(str));

%% dims is a cell array with allowed dimensions. 
% S is the size of the string dimension. S=0 when S does not exist
% dimA is a vector with the actual dimensions.
% dim is a string with the matching dimension
function [dim,S]=checkdim(Obj,dims,dimA)
dim=[]; S=0;
for jj=1:length(dims)
  dimS=dims{jj};
  if length(dimS)==1, dimS=[dimS 'I']; end; % 1D required, but Matlab is always 2D at least.
	dimR=getdim(Obj.API,dimS);
	if length(dimA)==length(dimR), % the same size?    
    if ~isempty(strfind(dimS,'S'))
      Sidx=strfind(dimS,'S');
      S=max(S,dimA(Sidx));
      dimR(Sidx)=dimA(Sidx); % string dim are always correct
    end
		if dimA==dimR, dim=upper(dims{jj}); break; end;	% found!
	elseif length(dimA)<length(dimR)	% extend the size?
		if [dimA ones(1,length(dimR)-length(dimA))]==dimR, dim=upper(dims{jj}); break; end; % found!
	end
end
