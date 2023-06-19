% this is a simple function to calculate tour length
function length=TourLength(tour,L)

    n=numel(tour); % length of the tour
    tour=[tour tour(1)]; % appending the start to the end to complete the tour
    length=0;
    for i=1:n
        length=length+L(tour(i),tour(i+1));
    end

end