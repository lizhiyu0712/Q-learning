function j=RouletteWheelSelection(P)

    % please see the comments to understand Roulette Wheel Selection
    % algorithm
    
    r=rand; % generating a random number b/w 0 to 1
    C=cumsum(P); % cumulative sum
    j=find(r<=C,1,'first'); % pick the index of the first non-zero element

end