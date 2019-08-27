%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference spectrum generator for Micro-focus x-ray tubes    %
%                                                             %
% INPUT PARAMETERS:                                           % 
%                                                             %
%    varargin{1}: Maximum energy of the spectrum. A valid     %
%    range between 10keV to 90keV.                            %
%                                                             %
%    varargin{2}: Normalize spectrum. A boolean value, 1      %
%    to enable normalization and 0 to disable normalization.  %
%                                                             %
%    varargin{3} ... varargin{n}: Add beam filteration        %
%    element symbol as a string followed by its thickness     %
%    in mm.                                                   %
%                                                             %
% EXAMPLE USAGE:                                              %
%                                                             %
%  [spec, trans] = MicrofocusSpec(90,1,'Al',0.5,'Yb',0.25)    %
%                                                             %
% OUTPUT:                                                     %
%                                                             %
%    varargout{1}: filtered or unfilterd spectrum depending   %
%    on the input filter arguments.                           %
%                                                             %
%    varargout{2}: filter transmission data.                  %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:                                                     %
%    Kunal Kumar,                                             %
%    Copyright, 2019                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright Notice: Permission to use, copy, modify, and distribute this 
% application and its components for educational, research, and not for 
% profit purposes, without a signed licensing agreement, is hereby granted,
% provided that the above copyright notice, this paragraph and the 
% following paragraph appear in all copies, modifications, and distributions.
%
% Disclaimer: This application has been developed solely for the purpose of
% research and education, use of this application for clinical dosimetry is
% strongly discouraged without a complete validation and regulatory 
% approval. In no event shall the author of this application be liable to 
% any party for direct, indirect, special, incidental, or consequential 
% damages arising out of the use of this application. 

function [varargout] = MicrofocusSpec(varargin)

disp("### info: valid energy range: 10keV to 90keV")

switch nargin
    case 0
        error('Please specify the maximum energy !!');
        
    case 1
        energy = varargin{1}; % keV
        spectrum = read_specData(energy, 0);
        varargout{1} =  spectrum;
        
    case 2
        energy = varargin{1}; % keV
        Normalize     = varargin{2}; % Boolean
        spectrum = read_specData(energy, Normalize);
        varargout{1} =  spectrum;
        
    otherwise
        if ischar(varargin{3}) == 0
            
            error('Third argument must be a string. Please enter a filter element !!');
            
        else
            energy = varargin{1}; % keV
            Normalize     = varargin{2}; % Boolean
            spectrum = read_specData(energy, 0);
            filters = varargin(3:2:nargin-1);
            thickness =  varargin(4:2:nargin);
            
            transmission = filter_Transmission( filters, thickness);
            
            trans = transmission(:,1);
            for i = 1:(size(transmission,2)-1)
                trans = transmission(:,i+1).*trans; % multiply columns
            end
            
            FilteredSpectra = spectrum.*trans;
            
            if Normalize == 1
            
                FilteredSpectra = NormalizeSpec(FilteredSpectra);
            end
            
            varargout{1} =  FilteredSpectra;
            varargout{2} =  transmission;
            
        end
end

end

function spec = read_specData(varargin)

energy = varargin{1};
Normalize =  varargin{2};

load 'tasmip.mat' tasmip;
spec = zeros(150,1);
for i=7:1:energy
    spec(i-1,1) = tasmip(i,1) + tasmip(i,2).*energy + tasmip(i,3).*(energy^2) + tasmip(i,4).*(energy^3);
end

if Normalize == 1

    spec = NormalizeSpec(spec);
end

end

function transmission = filter_Transmission( filters, thickness)

Energy = 1000:1000:150000;  % eV
t_nm = cell2mat(thickness).*1e+6;   % nm
MassAttCoeff = PhotonAttenuation(filters, Energy./10^6, 'mac'); % cm^2/g
Density = PhysProps(filters); % g/cm^3
LinAttCoeff = MassAttCoeff.*cell2mat(Density(:,2))'; % cm^-1
AttLength = 1./LinAttCoeff*10^7; % nm
transmission = exp(-t_nm./AttLength);

end

function NormalizedSpec = NormalizeSpec(spec)

  NormalizedSpec = spec./sum(spec);
end
