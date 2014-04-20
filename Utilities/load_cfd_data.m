function [CFDdata] = load_cfd_data(filename,gamma)
% Function for loading CFDdata and converting surface pressures to
% pressure coefficients if necessary
% This function is copied from NASA LaRC-Analytical Mechanics Associates
% supplied files as it contains additional data to be added to the MSL
% Cp look-up table
%
% Inputs
% filename = name of the file containing the CFD data
% gamma = specific heat ratio -- used to generate pressure
%
% Outputs
% CFDdata = structure that has Cp table from the CFD data


%Load CFDdata - Surface Pressures
disp(['Loading CFD Database: ' filename])
CFDdata = load(filename);

%CFD Freestream Conditions at each Mach point
%    [Mach# vinf(m/s) rhoinf(kg/m^3)]
if strcmp(filename(1:3),'MSL')
    cfd_addinfo =   [1.5,  327,   0.005273;
                    2.1,  472,  0.004634;
                    4.0,  892,  0.003834;
                    5.8,  1315, 0.003776;
                    8.8,  2000, 0.003483;
                    10.6, 2400, 0.003168;
                    16.1, 3600, 0.001926;
                    18.0, 4000, 0.0015166;
                    21.0, 4600, 0.0009331;
                    25.3, 5200, 0.0003448];
    
    %Convert to pressure coefficients
    % Pinf = (rhoinf/gamma)*(vinf/Mach)^2
    Pinf = (cfd_addinfo(:,3)/gamma).*(cfd_addinfo(:,2)./cfd_addinfo(:,1)).^2;
    % qbar = 0.5*rhoinf*Velinf^2
    qbar = .5*cfd_addinfo(:,3).*cfd_addinfo(:,2).^2;
    % Cp = (P - Pinf)/qbar
    for i = 1:size(cfd_addinfo,1)
        j = find(CFDdata.mach_vec == cfd_addinfo(i,1));
        CFDdata.mach(i).Cp = (CFDdata.mach(j).P-Pinf(i))/qbar(i);
    end
end