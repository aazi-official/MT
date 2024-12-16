function [AppRes,Phase] = generateSynData
% generate a synthetic MT data set
% Author: Jiajia Sun, 09/05/2015 at Colorado School of Mines
% Revision History:
% 12/14/2024: clean up the code for Xiaolong.
%==========================================================================

      res = [500,1000,100,500,1000];  % in ohm*meter
thickness = [210,1624,1346,1435];  % in meters

%        res = [500,500,500,500,500];  % in ohm*meter
% thickness = [210,1624,1346,1435];  % in meters

  NumFreq = 50;
frequency = logspace(3,-3,NumFreq);

AppRes = zeros(NumFreq,1);
Phase = zeros(NumFreq,1);

for j = 1:NumFreq
   [tmp1, tmp2,~] = MTmodeling1D(res,thickness,frequency(j));
   AppRes(j,1) = tmp1;
   Phase(j,1) = tmp2;
end

%Phase = Phase*180/pi;  % convert radians to degrees

% Add Gaussian noise
 % AppRes = log10(AppRes);
 % 
 % AppRes = AppRes + .01*randn(NumFreq,1);
 % Phase  = Phase  + .01*randn(NumFreq,1);
 
 %save -v7.3 AppRes.mat AppRes
 %save -v7.3 Phase.mat Phase

figure(1)
loglog(frequency,AppRes,'-ro');
set(gca,'XDir','reverse');
xlabel('frequencies');
ylabel('Apparent resistivities');

figure(2)
semilogx(frequency,Phase*180/pi,'-ro')
set(gca,'XDir','reverse');
xlabel('frequencies')
ylabel('Phases in degrees')


end





