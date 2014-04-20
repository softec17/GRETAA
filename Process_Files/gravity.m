function g = gravity(r,long,lat,mu,J2,gravmodel)
% Calculates the current gravitation acceleration
%
% Inputs
% r = planet-centric radius
% long = longitude (currently not used)
% lat = planet-centric latitude
% mu = planetary gravational constant
% J2 = planetary zonal harmonic of the second order constant
% gravmodel = flag for what kind of gravitional model
%
% Output
% g = gravitational acceleration vector in the vehicle-carried frame


switch gravmodel
    case 'rsquared'
        g = [0;0;mu/r^2];
%           g = [0;0;3.7*(3396.2e3)^2/r^2];
    case 'J2'
        g = [-3*mu*J2*sin(2*lat)/(2*r^4);...
            0;...
            mu/r^2-3*mu*J2*(2-3*cos(lat)^2)/(2*r^4)];
end

return