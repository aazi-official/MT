function [apparentRes, phase, impedances] = MTmodeling1D(res,thickness,frequency)
% 1D MT modeling program
%==========================================================================
% Input variables:
%       res: resistivity values for each layer from the toppest to the very
%            bottom layers
% thickness: the thickness of each layer. The very bottom layer is the
%            halfspace, therefore, its thickness is not given. That is, the
%            dimension of the vector thickness is one less than that of the
%            vector res.
% frequency: the frequency that is under study
%==========================================================================
% Output variables:
% apparentRes: Apparent resistivity at the given frequency
%       phase: the phase in radians (not in degrees)
%  impedances: impedance at each layer
%==========================================================================
% Author: Jiajia Sun, Sept 4th, 2015 at Colorado School of Mines
% Revision history:
%
%
%==========================================================================

mu = 4*pi*1E-7;   % Magnetic permeability (H/m)
omega = 2*pi*frequency;  % Angular frequency
NumLayer = numel(res);   % Number of layers

impedances = zeros(NumLayer,1);  % impedance
Zn = sqrt(sqrt(-1)*omega*mu*res(NumLayer)); % impedance of the halfspace
impedances(NumLayer,1) = Zn;

for j = NumLayer-1:-1:1
    Alphaj = sqrt(sqrt(-1)*omega*mu/res(j)); % inducation parameter
    Wj = sqrt(sqrt(-1)*omega*mu*res(j)); % intrinsic impedance
    Expj = exp(-2*Alphaj*thickness(j));  % Exponential factor
    
    BelowImpedance = impedances(j+1);
    Rj = (Wj - BelowImpedance)/(Wj + BelowImpedance); % reflection coefficient
    REj = Expj*Rj;
    Zj = Wj*((1-REj)/(1+REj)); % impedance
    impedances(j,1) = Zj;    
end

Z = impedances(1);
apparentRes = abs(Z)^2/(omega*mu);
phase = atan(imag(Z)/real(Z));

end




















