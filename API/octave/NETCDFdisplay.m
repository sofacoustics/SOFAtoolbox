function [] = NETCDFdisplay(filename)
%NETCDFDISPLAY 
%   [] = NETCDFdisplay(filename) displays information about specified SOFA file

% SOFA API - function octave/NETCDFdisplay
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

ncid=netcdf(filename,'r');
temp=ncvar(ncid);
for ii=1:length(temp)
    fprintf(['--- ' ncname(temp{ii}) ' ---\n'])
    fprintf(['DataType: ' ncdatatype(temp{ii}) '\n'])
    fprintf('Dimensions:\n')
    for jj=1:length(ncdim(temp{ii})(:))
        disp(ncdim(temp{ii}){jj})
    end
    fprintf('\n')
end

end %of function