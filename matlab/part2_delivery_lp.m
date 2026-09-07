% Recovered Part 2: transportation linear program.
%
% The original part2.mat is missing, so this script reads the editable CSV
% reconstruction in data/provinces.csv.

clear; clc;

T = readtable(fullfile('data', 'provinces.csv'), 'TextType', 'string');
n = height(T);

p = T.population_index;
e = T.economic_index;
b = T.supply_capacity;
c_demand = T.demand;

if sum(b) < sum(c_demand)
    b = b * (sum(c_demand) / sum(b)) * 1.02;
end

D = distance_index_from_lat_lon(T.latitude, T.longitude, 6, 4);

cost = zeros(n, n);
for i = 1:n
    for j = 1:n
        cost(i, j) = p(j) * ((e(i) + e(j)) / 2) * D(i, j);
    end
end

f = reshape(cost', [], 1);

% Row sums <= supply.
A = zeros(n, n * n);
for i = 1:n
    for j = 1:n
        A(i, sub2ind([n, n], j, i)) = 1;
    end
end

% Column sums >= demand, represented as -column_sum <= -demand.
A_demand = zeros(n, n * n);
for j = 1:n
    for i = 1:n
        A_demand(j, sub2ind([n, n], j, i)) = -1;
    end
end

A = [A; A_demand];
rhs = [b; -c_demand];

Aeq = [];
beq = [];
lb = zeros(n * n, 1);

options = optimoptions('linprog', 'Display', 'off');
[weight, total_cost, exitflag] = linprog(f, A, rhs, Aeq, beq, lb, [], options);

if exitflag <= 0
    error('linprog did not find an optimal solution.');
end

X = reshape(weight, n, n)';

fprintf('Reconstructed Part 2 cost: %.4f\n', total_cost);
disp(array2table(X, 'VariableNames', cellstr(T.capital), 'RowNames', cellstr(T.capital)));

function D = distance_index_from_lat_lon(lat, lon, low, width)
    n = numel(lat);
    km = zeros(n, n);
    for i = 1:n
        for j = 1:n
            km(i, j) = haversine_km(lat(i), lon(i), lat(j), lon(j));
        end
    end
    D = low + width * km / max(km, [], 'all');
end

function d = haversine_km(lat1, lon1, lat2, lon2)
    r = 6371;
    phi1 = deg2rad(lat1);
    phi2 = deg2rad(lat2);
    dphi = deg2rad(lat2 - lat1);
    dlambda = deg2rad(lon2 - lon1);
    h = sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlambda / 2)^2;
    d = 2 * r * asin(sqrt(h));
end
