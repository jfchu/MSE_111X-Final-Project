%% PROBLEM 1
% Assign values 
k = 50;
m = 10;
n = 10000;

b = 1000*ones(m, 1); 

% Offline problem
% To allow for meaningful comparison, Problem 1-3 will draw from the same
% randomly generated constraints matrix and bid vector  

% m-by-n matrix of 0s and 1s
A = randi([0 1], m, n);

% ground truth price vector whatever that is
p = ones(m, 1);

% sequence of random bids 
bids = [];

for j = 1:n
    % randn(0, 0.2) is a lie, only accepts integer input. Damn
    newbid = dot(p, A(:, j)) + 0.2*randn;
    bids = [bids; newbid];
end

% Find primal solution for offline case
cvx_begin
    variables x(n);
    maximize dot(bids, x);
    subject to
        A*x <= b;
        0 <= x <= 1;
cvx_end 
offline_opt = cvx_optval;

% Store solutions
sol_offline = x;

% Online problem 
% Index bids and constraints from offline problemn

A_online = A(:, 1:k);
bids_online = bids(1:k);

% Find primal and dual optimal solutions for online case
cvx_begin
    variables x(k);
    dual variable y_1;
    dual variable y_2;
    dual variable y_3;
    maximize dot(bids_online, x);
    subject to
        y_1: A_online*x <= (k/n)*b;
        y_2: x >= 0;
        y_3: x <= 1;
cvx_end 

% Store optimal solution
solution = x;

% We use 'solution' later, prevent override
sol_50 = solution;

% Fixed dual price
dual_50 = y_1;

% Siumlate and store revenue from k + 1 to n bids
% We do not subtract from b yet, so always remaining goods left
for i = k + 1:n
    if dot(A(:, i), dual_50) < bids(i)
        sol_50 = [sol_50; 1];
    else 
        sol_50 = [sol_50; 0];
    end
end  

online_opt = dot(bids, sol_50);

% Find and store primal and dual optimal solutions for k = 100
cvx_begin
    variables x(100);
    dual variable y_1;
    dual variable y_2;
    dual variable y_3;
    maximize dot(bids(1:100), x);
    subject to
        y_1: A(:, 1:100)*x <= (100/n)*b;
        y_2: x >= 0;
        y_3: x <= 1;
cvx_end 
sol_100 = x;
dual_100 = y_1;

% Siumlate and store revenue from 101 to n bids
for i = 101:n
    if dot(A(:, i), dual_100) < bids(i)
        sol_100 = [sol_100; 1];
    else 
        sol_100 = [sol_100; 0];
    end
end  

online_opt_2 = dot(bids, sol_100);

% Find and store primal and dual optimal solutions for k = 200
cvx_begin
    variables x(200);
    dual variable y_1;
    dual variable y_2;
    dual variable y_3;
    maximize dot(bids(1:200), x);
    subject to
        y_1: A(:, 1:200)*x <= (200/n)*b;
        y_2: x >= 0;
        y_3: x <= 1;
cvx_end 
sol_200 = x;
dual_200 = y_1;

% Siumlate and store revenue from 201 to n bids
for i = 201:n
    if dot(A(:, i), dual_200) < bids(i)
        sol_200 = [sol_200; 1];
    else 
        sol_200 = [sol_200; 0];
    end
end  

online_opt_3 = dot(bids, sol_200);

% Take quotients between online and offline optimal value 
value_quotient = online_opt/offline_opt;
value_quotient_2 = online_opt_2/offline_opt;
value_quotient_3 = online_opt_3/offline_opt;

% FIRST PLOT: quotients between online and offline optimal values 

% x-axis: k = 50, 100, 200
ratio_label = categorical({'50', '100', '200'});
ratio_label = reordercats(ratio_label,{'50', '100', '200'});

% y-axis: the quotients 
ratio = [value_quotient; value_quotient_2; value_quotient_3];

figure(1)
bar(ratio_label, ratio);
title('Ratios Between Online and Offline Optimal Values', 'interpreter', 'latex');
xlabel('first $k$ bids', 'interpreter', 'latex')
ylabel('$\bar{y}^k/\bar{y}^n$', 'interpreter', 'latex')

% SECPMD PLOT: accumulation of offline profit across bids
% Primairly for demonstrating linearity
figure(2)
plot(length, off_solution_sum);
title('Accumulation of Offline Algorithm Profit Across Bids', 'interpreter', 'latex');
xlabel('number of bids', 'interpreter', 'latex');
ylabel('current offline profit', 'interpreter', 'latex');

% THIRD PLOT: difference between accumulation of profit across bids,
% online solution (k = 50, 100, 200) and offline solution 
% Acts as further commentary on quotients graph 

% x-axis: 1 to n
length = [];
for i = 1:n
    length = [length; i];
end

% accumulation of offline profit
off_solution_sum = [];
for i = 1:n
    sum_append = dot(bids(1:i), sol_offline(1:i));
    off_solution_sum = [off_solution_sum; sum_append];
end

% accumulation of online profit (k = 50)
solution_sum_50 = [];
for i = 1:n
    sum_append = dot(bids(1:i), sol_50(1:i));
    solution_sum_50 = [solution_sum_50; sum_append];
end

% accumulation of online profit (k = 100)
solution_sum_100 = [];
for i = 1:n
    sum_append = dot(bids(1:i), sol_100(1:i));
    solution_sum_100 = [solution_sum_100; sum_append];
end

% accumulation of online profit (k = 200)
solution_sum_200 = [];
for i = 1:n
    sum_append = dot(bids(1:i), sol_200(1:i));
    solution_sum_200 = [solution_sum_200; sum_append];
end

figure(3)
plot(length, solution_sum_50 - off_solution_sum);

hold on
plot(length, solution_sum_100 - off_solution_sum);

hold on
plot(length, solution_sum_200 - off_solution_sum);

title('Difference Between Online and Offline Algorithm Profit Across Bids', 'interpreter', 'latex');
xlabel('number of bids', 'interpreter', 'latex');
ylabel('$\bar{y}^k$ minus offline profit', 'interpreter', 'latex');
legend('$k = 50$', '$k = 100$', '$k = 200$', 'interpreter', 'latex');

%% PROBLEM 2
% Initiate matrix to keep track of dual solutions
B = [];

% Always set k to 50 for this problem 
while k <= n
    
    % Find dual optimal solution at k bids 
    cvx_begin
        variables x(k);
        dual variable y_1;
        dual variable y_2;
        dual variable y_3;
        maximize dot(bids_online, x);
        subject to
            y_1: A_online*x <= (k/n)*b;
            y_2: x >= 0;
            y_3: x <= 1;
    cvx_end 

    % Build matrix of duals 
    B = [B y_1];

    % Update dual 
    dual_part2 = y_1;
    
    k = 2*k;

    % Iteration starts from 51, 101, . . . 3201, 6401
    first_iterate = k/2 + 1;

    % Iteration ends at 100, 200, . . . 6400, 10000
    last_iterate = min([k, n]);

    % Update indexed bids and constraints from offline problem
    A_online = A(:, 1:last_iterate);
    bids_online = bids(1:last_iterate);
   
    % Iterate decisions based on updated dual until certain time points
    for i = first_iterate:last_iterate
        if dot(A_online(:, i), dual_part2) < bids_online(i)
            solution = [solution; 1];
        else 
            solution = [solution; 0];
        end
    end       
end

dynamic_opt = dot(bids, solution);

% FOURTH PLOT: accumulation of dynamic online revenue across bids
% Primairly for demonstrating linearity 
figure(4);
plot(length, solution_sum);
title('Accumulation of Dynamic Online Algorithm Profit Across Bids', 'interpreter', 'latex');
xlabel('number of bids', 'interpreter', 'latex');
ylabel('current dynamic online profit', 'interpreter', 'latex');


% FIFTH PLOT: difference between accumulation of profit across bids,
% dynamic online and static online/offline algorithms

% accumulation of dynamic online solution
solution_sum = [];
for i = 1:n
    sum_append = dot(bids(1:i), solution(1:i));
    solution_sum = [solution_sum; sum_append];
end

figure(5)
plot(length, solution_sum - off_solution_sum);

hold on
plot(length, solution_sum - solution_sum_50);

hold on 
plot(length, solution_sum - solution_sum_100);

hold on
plot(length, solution_sum - solution_sum_200);

title('Difference Between Dynamic and Static Algorithm Profit Across Bids', 'interpreter', 'latex');
xlabel('number of bids', 'interpreter', 'latex');
ylabel('dynamic online profit minus . . . ', 'interpreter', 'latex');
legend('offline', '$k = 50$', '$k = 100$', '$k = 200$', 'interpreter', 'latex')


% SIXTH PLOT: convergence of dual solution over time points
% Augment offline dual price for comparison
B = [B dual_solution];

% Distance between truth vector and dual prices over time points
N = [];
c = log2(k/50) + 1;
for i = 1:c
    g = p - B(:, i);
    N = [N; norm(g)];
end

% x-axis: k = 50, 100, 200, . . . 
time_points = categorical({'50', '100', '200', '400', '800', '1600', '3200', '6400', '10000'});
time_points = reordercats(time_points,{'50', '100', '200', '400', '800', '1600', '3200', '6400', '10000'});

figure(6)
bar(time_points, N);
title('Convergence of Dual Solution Across Time Points', 'interpreter', 'latex')
xlabel('number of bids', 'interpreter', 'latex');
ylabel('$\sqrt{(\bar{y}^k - \bar{p})^2}$', 'interpreter', 'latex');

%% PROBLEM 3

A_part3 = A(:, 1);

% Find optimal dual solution for single bid - initial case
cvx_begin
    variable x(1);
    dual variable y_1;
    dual variable y_2;
    dual variable y_3;
    maximize bids(1)*x;
    subject to
        y_1: A_part3*x <= 1/(n - 1)*b;
        y_2: x >= 0;
        y_3: x <= 1;
cvx_end

% Store dual optimal solution
dual_part3 = y_1;

% Initiate solution to action-history-dependent learning algorithm (AHD)
LP_decisions = [];

% Initiate vector of differences between profit at k iterations and and k/n
% times offline optimal value, or offline performance vector 
profit_off_diff = [];

% Initiate vector of differences between profit at k iterations and and k/n
% times dynamic online optimal value, or online performance vector
profit_on_diff = [];

% Initiate vector of sum of b over number of bids 
% Acts as visualization of depletion of resources 
remaining_resources = [sum(b)];

if dot(A_part3(:, 1), dual_part3) < bids(1)
    LP_decisions = [LP_decisions; 1];
else
    LP_decisions = [LP_decisions; 0];
end

% s is equivalent to k for previous problems, but don't want to risk
% overriding results from previous parts 

% Keep track of iterations
i = 1;

% AHD start
% iteration stops right before s = n to avoid dividing by 0
for s = 2:n-1

    % Index relavant part of matrix
    A_part3 = A(:, 1:s);
    
    % Augment based on dual at sth bid
    % min(b) - if one element within b is negative, we've already run out
    % of resources
    if dot(A_part3(:, s), dual_part3) < bids(s) && min(b) >= 0
        LP_decisions = [LP_decisions; 1];
    else 
        LP_decisions = [LP_decisions; 0];
    end
    
    i = i + 1;

    % Increment offline performance vector
    profit_off_diff_append = dot(bids(1:s), LP_decisions) - s/n*offline_opt;
    profit_off_diff = [profit_off_diff; profit_off_diff_append];

    % Increment online performance vector
    profit_on_diff_append = dot(bids(1:s), LP_decisions) - s/n*dynamic_opt;
    profit_on_diff = [profit_on_diff; profit_on_diff_append];

    % Update number of resources
    b = b - LP_decisions(s)*A_part3(:, s);
   
    % Store sum of entries in updated resources vector 
    remaining_resources = [remaining_resources; sum(b)];
    
    % Break after updating necessary vectors to keep track of number of
    % iterations easier
    if min(b) < 0
        break;
    end

    % Update optimal dual solution for s bids
    cvx_begin
        variables x(s);
        dual variable y_1;
        dual variable y_2;
        dual variable y_3;
        maximize dot(bids(1:s), x);
        subject to
            y_1: A_part3*x <= s/(n - s)*b;
            y_2: x >= 0;
            y_3: x <= 1;
    cvx_end

    % Update dual prices
    dual_part3 = y_1;
end

% SEVENTH PLOT: performance of offline algorithm in AHD

figure(7)
t = tiledlayout(1, 2);
title(t, 'Offline vs. Dynamic Online Performance in AHD Learning Algorithm', 'interpreter', 'latex');
xlabel(t, 'number of bids', 'interpreter', 'latex');
ylabel(t, '$\sum_{j = 1}^{k} \pi_jx_j - \frac{k}{n}OPT$', 'interpreter', 'latex');

nexttile
plot(length(2:i), profit_off_diff);
xlabel('$OPT = \mbox{offline optimal value}$', 'interpreter', 'latex');

nexttile
plot(length(2:i), profit_on_diff);
xlabel('$OPT = \mbox{dynamic optimal value}$', 'interpreter', 'latex');
        
