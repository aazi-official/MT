data=load("output_langevin2_3.mat");
data=data.output_langevin1.samples;
%mean2=mean(output_tmcmc.samples);
mean2=mean(data);
%var2=var(output_tmcmc.samples);
var2=var(data);
result=[500, 1000, 100, 500, 1000, 200,210, 1624, 1346, 1435, 1800];  


rh = [500, 1000, 100, 500, 1000, 200];
thickness = [-210, -1624, -1346, -1435, -1800, -Inf]; 
std_dev2 = sqrt(var2); 


rh_result = result(1:6); 
thickness_result = [-result(7:end), -10000]; 

rh_mean2 = mean2(1:6); 
thickness_mean2 = [-mean2(7:end), -10000]; 

rh_upper = rh_mean2 + std_dev2(1:6); 
rh_lower = rh_mean2 - std_dev2(1:6); 

thickness_upper = [-mean2(7:end) + std_dev2(7:end), -1000000]; 
thickness_lower = [-mean2(7:end) - std_dev2(7:end), -1000000];

% stair
% Result 
thickness_stairs_result = repelem(thickness_result, 2);
thickness_stairs_result = thickness_stairs_result(1:end-1);
thickness_stairs_result = [0, thickness_stairs_result];
for i = 2:5
    thickness_stairs_result(2*i) = thickness_stairs_result(2*i-1) + thickness_stairs_result(2*i);
    thickness_stairs_result(2*i+1) = thickness_stairs_result(2*i);
end
rh_stairs_result = repelem(rh_result, 2);

% Mean2
thickness_stairs_mean2 = repelem(thickness_mean2, 2);
thickness_stairs_mean2 = thickness_stairs_mean2(1:end-1);
thickness_stairs_mean2 = [0, thickness_stairs_mean2];
for i = 2:5
    thickness_stairs_mean2(2*i) = thickness_stairs_mean2(2*i-1) + thickness_stairs_mean2(2*i);
    thickness_stairs_mean2(2*i+1) = thickness_stairs_mean2(2*i);
end
rh_stairs_mean2 = repelem(rh_mean2, 2);

% bound
thickness_stairs_upper = repelem(thickness_upper, 2);
thickness_stairs_upper = thickness_stairs_upper(1:end-1);
thickness_stairs_upper = [0, thickness_stairs_upper];
for i = 2:5
    thickness_stairs_upper(2*i) = thickness_stairs_upper(2*i-1) + thickness_stairs_upper(2*i);
    thickness_stairs_upper(2*i+1) = thickness_stairs_upper(2*i);
end
thickness_stairs_lower = repelem(thickness_lower, 2);
thickness_stairs_lower = thickness_stairs_lower(1:end-1);
thickness_stairs_lower = [0, thickness_stairs_lower];
for i = 2:5
    thickness_stairs_lower(2*i) = thickness_stairs_lower(2*i-1) + thickness_stairs_lower(2*i);
    thickness_stairs_lower(2*i+1) = thickness_stairs_lower(2*i);
end
rh_stairs_upper = repelem(rh_upper, 2);
rh_stairs_lower = repelem(rh_lower, 2);


figure;
stairs(rh_stairs_result, thickness_stairs_result, 'b-', 'LineWidth', 1.5, 'DisplayName', 'True Values'); 
hold on;
stairs(rh_stairs_mean2, thickness_stairs_mean2, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Predicted Mean');


% stairs(rh_stairs_upper, thickness_stairs_upper, 'r--', 'LineWidth', 1, 'DisplayName', 'Upper Bound (Predicted)');
% stairs(rh_stairs_lower, thickness_stairs_lower, 'c--', 'LineWidth', 1, 'DisplayName', 'Lower Bound (Predicted)');
x_fill = [rh_stairs_upper, fliplr(rh_stairs_lower)];
y_fill = [thickness_stairs_upper, fliplr(thickness_stairs_lower)];
fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Uncertainty Region');

xlabel('Resistivity (\Omega\cdotm)');
ylabel('Thickness (m)');
title('Comparison of True Values and LTMCMC Predicted Values with Uncertainty');
grid on;
xlim([0, 1200]);
ylim([-7000, 500]);
legend('Location', 'best');
hold off;