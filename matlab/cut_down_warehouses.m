% Recovered "Code for Cutting Down Warehouses".
%
% This follows the algorithm visible in the original Word-document screenshots:
% remove one warehouse at a time, redistribute its supply equally among the
% remaining active warehouses, solve the LP, and keep the cheapest removal.

clear; clc;

T = readtable(fullfile('data', 'provinces.csv'), 'TextType', 'string');
n = height(T);

p = T.population_index;
e = T.economic_index;
b_old = T.supply_capacity;
c_demand = T.demand;

if sum(b_old) < sum(c_demand)
    b_old = b_old * (sum(c_demand) / sum(b_old)) * 1.02;
end

D = distance_index_from_lat_lon(T.latitude, T.longitude, 6, 4);
base_cost = zeros(n, n);
for i = 1:n
    for j = 1:n
        base_cost(i, j) = p(j) * ((e(i) + e(j)) / 2) * D(i, j);
    end
end
f = reshape(base_cost', [], 1);

[A, rhs_demand, lb] = build_constraints(n, c_demand);

removed_indices = [];
removed_costs = [];
target_warehouses = 16;

while sum(b_old > 0) > target_warehouses
    total_costs = inf(n, 1);
    candidate_supplies = cell(n, 1);

    active = find(b_old > 0);
    for idx = 1:numel(active)
        i = active(idx);
        b = b_old;
        receivers = active(active ~= i);
        b_temp = b(i) / numel(receivers);
        b(i) = 0;
        b(receivers) = b(receivers) + b_temp;

        rhs = [b; rhs_demand];
        options = optimoptions('linprog', 'Display', 'off');
        [~, total_cost, exitflag] = linprog(f, A, rhs, [], [], lb, [], options);

        if exitflag > 0
            total_costs(i) = total_cost;
            candidate_supplies{i} = b;
        end
    end

    [m, k] = min(total_costs);
    if isinf(m)
        error('No feasible warehouse removal found.');
    end

    b_old = candidate_supplies{k};
    removed_indices = [removed_indices; k]; %#ok<AGROW>
    removed_costs = [removed_costs; m]; %#ok<AGROW>

    fprintf('Removed %s (%s), active warehouses: %d, cost: %.4f\n', ...
        T.province(k), T.capital(k), sum(b_old > 0), m);
end

fprintf('\nFinal active capitals:\n');
disp(T.capital(b_old > 0));

history = table((1:numel(removed_indices))', removed_indices, ...
    T.province(removed_indices), T.capital(removed_indices), removed_costs, ...
    'VariableNames', {'step', 'removed_index_1_based', 'removed_province', 'removed_capital', 'total_cost'});
disp(history);

function [A, rhs_demand, lb] = build_constraints(n, demand)
    A_supply = zeros(n, n * n);
    for i = 1:n
        for j = 1:n
            A_supply(i, sub2ind([n, n], j, i)) = 1;
        end
    end

    A_demand = zeros(n, n * n);
    for j = 1:n
        for i = 1:n
            A_demand(j, sub2ind([n, n], j, i)) = -1;
        end
    end

    A = [A_supply; A_demand];
    rhs_demand = -demand;
    lb = zeros(n * n, 1);
end

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
