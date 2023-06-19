
%% -------- STEP 0 --------

clc;
clear;
close all;

% Our goal is to find the optimal path amongst few locations, while 
% visiting each node only once. we will use ant colony optimization
% algorithm to solve this.

%% -------- STEP 1 --------

% setup the problem
% first we will set the location where the salesperson needs to visit, in a
% X-Y coordinate system.

% n = 16; % Total number of locations  
% x = round(100*rand(1,n));
% y = round(100*rand(1,n));

n = 16; % total number of locations  
% Coordinates of the locations
x = [36    13    76    99    27    62    94    74    31    67    56    39    66    44    20    55];
y = [59    76    38    89    66    83     5    50    42    61    80    10    87    67    64    32];


% Visualize the locations
figure(1)
scatter(x,y, 100, 'MarkerFaceColor','y', 'MarkerEdgeColor','k','LineWidth',1.5)
xlabel('x');
ylabel('y');
axis equal;
grid on;

% Initialize the distance matrix
L = zeros(n, n);
% Fillup the distance matrix
for i = 1:n-1
    for j = i+1:n
        L(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
        L(j,i) = L(i,j);
    end
end

%% -------- STEP 2 --------

% Ant Colony Optimization (ACO) setup

% Parameters
max_iter = 1;  % Maximum number of iteration
num_ants = 10;    % Number of ants
Q        = 1;     % q constant
alpha    = 1;     % Phromone coefficient
beta     = 1;     % Preference for shortest path
rho      = 0.1;   % Evaporation rate
tau0     = 10*Q/(n*mean(L(:)));	% Initial phromone


% Initilize parameters
eta           = 1./L; % Desirability
tau           = tau0*ones(n,n); % Phromone matrix
best_cost     = zeros(max_iter,1); % We are creating an array to hold the best cost values
best_sol.cost = inf; % Best solution holder

% Ant colony matrix
empty_ant.tour = [];
empty_ant.cost = [];
ant = repmat(empty_ant,num_ants,1);

% Performing the search operation
for it = 1:max_iter % Stop when maximum iteration reached 
    for k = 1:num_ants % Iterate over all the ants
        
        % starting at a random location
        ant(k).tour = randi([1 n]);
        
        for l = 2:n % iterate of all the locations
            
            % Last location of the travel path
            i = ant(k).tour(end);
            % Calculating probability P_xy for kth ant; Slide 95 of the lectures
            P = tau(i,:).^alpha.*eta(i,:).^beta;
            P(ant(k).tour) = 0;
            P = P/sum(P); 
            
            % We will use the calculated probabilities in 
            % Roulette Wheel Selection method to select the next location
            
            j = RouletteWheelSelection(P); % the next location
            ant(k).tour = [ant(k).tour j];
            
        end
        
        % Updating the cost function for kth ant
        ant(k).cost = TourLength(ant(k).tour, L);
        
        if ant(k).cost < best_sol.cost
            best_sol = ant(k);  % Updating the best solution so far
        end
        
    end
    
    % Update Phromones
    for k = 1:num_ants
        
        tour = ant(k).tour; % Travel path for ant k
        tour = [tour tour(1)]; % Appending the start to the end to complete the tour
        
        for l=1:n
            
            i = tour(l);
            j = tour(l+1);
            
            % Updating the phromones with the cost following the equation on Slide 90
%             tau(i,j) = (1-rho)*tau(i,j)+Q/ant(k).cost;
            
        end
    end
    
    % update best cost
    best_cost(it) = best_sol.cost;
    
    % printing cost updates
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(best_cost(it))]);
        
end

% Results
final_tour = [best_sol.tour best_sol.tour(1)];
disp(['best path: ' num2str(final_tour)])


%% Visualizing the travel path

figure(2)
plot(x(final_tour),y(final_tour),'k-o',...
        'MarkerSize',10,...
        'MarkerFaceColor','y',...
        'LineWidth',1.5);
xlabel('x');
ylabel('y');
axis equal;
grid on;
    
% Visualize the cost curve
figure;
plot(best_cost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;




