% Recovered Part 1: minimize the number of additional warehouses.
%
% This is a cleaned transcription of the MATLAB code screenshot in the paper.
% It uses the 14 "non" provinces. Recovery note: when solved literally,
% the screenshot code does not reproduce the six selected locations reported
% in the paper text. Both outputs are printed below.

clear; clc;

candidate_provinces = { ...
    'Inner Mongolia', 'Xinjiang', 'Tibet', 'Ningxia', ...
    'Heilongjiang', 'Liaoning', 'Hebei', 'Henan', ...
    'Gansu', 'Shaanxi', 'Qinghai', 'Zhejiang', ...
    'Jiangxi', 'Yunnan'};

f = ones(14, 1);
intcon = 1:14; % every index is an integer

% Original screenshot stores the constraint as A*x <= b with negative values.
% This is equivalent to rent_values*x >= 386.9.
b = -386.9 * ones(14, 1);
A_row = [-31.66, -29.5, -31.99, -30.48, -23.95, ...
         -22.04, -22.77, -23.43, -17.93, -20.74, ...
         -19.63, -21.07, -19.2, -150];
A = repmat(A_row, 14, 1);

Aeq = [];
beq = [];
lb = zeros(14, 1);
ub = ones(14, 1);

options = optimoptions('intlinprog', 'Display', 'off');
[x, objective_value, exitflag] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, options);

if exitflag <= 0
    error('intlinprog did not find an optimal solution.');
end

selected = candidate_provinces(logical(round(x)));

fprintf('Minimum extra warehouses: %.0f\n', objective_value);
fprintf('Selected provinces from literal recovered code:\n');
disp(selected');

fprintf('Paper reported selected provinces:\n');
disp({'Zhejiang'; 'Liaoning'; 'Henan'; 'Heilongjiang'; 'Jiangxi'; 'Gansu'});
