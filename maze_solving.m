
%% -------- STEP 0 --------

clc;
clear;close all;



%% -------- STEP 1 --------

% Maze solving using Q-learning algorithm

% First create a maze. When doing so, please make sure it is a
% solvable one, i.e., if the start or finish states were blocked, 
% you will need to re-run the program.

% Creating an n*n maze 
n = 12;

% obst is the number of obstacles in the maze
obst = 30;

maze = ones(n,n);
for i = 1:obst
    maze(randi([1,n]),randi([1,n])) = -100;
end
maze(1,1) = 1; % start point
maze(n,n) = 10; % finish point


% Visualize the maze in matrix form
disp(maze)

% Visualize the maze as a figure
n = length(maze);
figure
imagesc(maze)
mycolors = [0.4 0.4 0.4; 0.8 0.8 0.8];
colormap(mycolors);

for i=1:n
    for j = 1:n
        if maze(i,j) == -100
            text(j,i,'X','HorizontalAlignment','center')
        end
    end
end
text(1,1,'START','HorizontalAlignment','center')
text(n,n,'FINISH','HorizontalAlignment','center')
axis off

START = 1;
FINISH = n*n;
% fprintf('Final State is: %d\n',FINISH)


%% -------- STEP 2 --------

% Next we will create the reward matrix for the maze
% First, let's desfine the possible actions that the agent can take in this
% environment. 
% 
% Let's assume at each state we have 8 possible actions:
% upward                    :  (i-n)
% downward                  :  (i+n)
% left                      :  (i-1)
% right                     :  (i+1)
% diagonally south east     :  (i+n+1)
% diagonally south west     :  (i+n-1)
% diagonally north east     :  (i-n+1)
% diagonally north west     :  (i-n-1)
% For any other action it should receive reward of -INF (negative infinity) 
% to prevent it from happening

reward = zeros(n*n);
for i = 1:FINISH
    reward(i,:) = reshape(maze',1,FINISH);
end
for i = 1:FINISH
    for j = 1:FINISH
        % This line enables the permitted actions
        if j ~= i-n  && j ~= i+n  && j ~= i-1 && j ~= i+1 
            reward(i,j)=-Inf;
        end    
    end
end
for i = 1:n:FINISH
    for j = 1:i+n
        if j == i+n-1 || j == i-1 || j == i-n-1
            reward(i,j) = -Inf;
            reward(j,i) = -Inf;
        end
    end
end


%% -------- STEP 3 --------

% Next, we aim to perform Q-learning.

% Recall from your lecture notes that we have a number of parameters that
% need to be set, namely: the Q-matrix, discount factor (gamma), learing
% rate (alpha), maximum number of iteration (max_iter).

qtable = randn(size(reward)); % Initialize randomly
gamma = 0.99; 
max_iter = 1000;
alpha = 0.1;

% We will find the optimal policy using the Bellman's equation; 
for i = 1:max_iter % Stop when maximum number of iterations reached
    
    % cs means current state
    cs = START;
    while(1) % Repeat untill the goal is reached

        % Possible actions for the chosen state
        n_actions = find(reward(cs,:)>0);

        % Choose a random action and set it as the next state
        % ns means next state
        ns = n_actions(randi([1 length(n_actions)],1,1)); 
       
        % Find all the possible actions for the selected state
        n_actions = find(reward(ns,:)>=0);

        % Find the maximum q-value i.e., the next state with the best action
        max_q = 0;
        for j = 1:length(n_actions)
            max_q = max(max_q,qtable(ns,n_actions(j)));
        end

        % Here we update the q-values as per Bellman's equation
        qtable(cs,ns) = (1-alpha) * qtable(cs,ns) + alpha * (reward(cs,ns) + gamma * max_q);

        % Check whether the episode has completed (reached the FINISH point)
        if(cs == FINISH)
            break;
        end

        % Set current state as next state
        cs = ns;
    end
end


%% -------- STEP 4 --------

% Finally we will solve and visualize the maze

start = START;
move = 0;
path = start;

% Iterate till the FINISH point is reached
while(move ~= FINISH)
    [~,move] = max(qtable(start,:));
    
    % Deleting chances of getting stuck in small loops (upto order of 4)  
    if ismember(move,path)
        [~,x] = sort(qtable(start,:),'descend');
        move = x(2); 
        if ismember(move,path)
            [~,x] = sort(qtable(start,:),'descend');
            move = x(3);
            if ismember(move,path)
                [~,x] = sort(qtable(start,:),'descend');
                move = x(4);
                if ismember(move,path)
                    [~,x] = sort(qtable(start,:),'descend');
                    move = x(5);
                end
            end
        end
    end
    
    % Appending next action/move to the path
    path  = [path,move];
    start = move;
end


% Reproducing the path to the matrix path
pmat    = zeros(n,n);
[q, r]  = quorem(sym(path),sym(n));
q       = double(q);
r       = double(r);
q(r~=0) = q(r~=0)+1;
r(r==0) = n;

for i = 1:length(q)
    pmat(q(i),r(i)) = 100;
end  


% Result: optimal path between START to FINISH
fprintf('Total steps: %d\n',length(path))
fprintf('Final path:\n')
for row = 1:length(pmat)
    for col = 1:length(pmat)
        if pmat(row,col) == 100
            fprintf('(%d,%d)\t', row,col);
        end
    end
end


% Visualize the solution
figure
imagesc(pmat)

p = 1;
for i = 1:n
    for j = 1:n
        if maze(i,j) == -100
            text(j,i,'X','HorizontalAlignment','center')
        end
        if pmat(i,j) == 100
            text(j,i,'\bullet','Color','black','FontSize',28)
            p=p+1;
        end
    end
end

text(1,1,'START','HorizontalAlignment','right')
text(n,n,'GOAL','HorizontalAlignment','right')

hold on
imagesc(maze,'AlphaData',0.2)
mycolors = [0.1 0.1 0.1; 0.8 0.8 0.8; 0.8 0.4 0.8];
colormap(mycolors);

hold off
axis off



