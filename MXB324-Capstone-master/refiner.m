function d = refiner(N, L, refiner_weights, refinements)
%REFINER Makes a refined one-dimensional mesh for the given N, L, and
%refinement ranges.
%
%Outputs: 
%- d = The deltas along this axis as a function of node index on this axis.
%
%Inputs:
%- N = The number of nodes that we will have along this axis.
%- L = The length of this axis in metres.
%- refiner_weights = A pair of values representing the smallest refinement
%  weight and the largest refinement rate respectively. Make sure 
%- refinements = An array of each refined point, which are themselves
%  a pair of values (the beginning and ending index of refined area) and
%  then a third element which is either 0, 1, or 2 which specifies if this
%  is at a min, mid, or max of the axis.
%
%Task name: mesh
%Authors: Paul
%
%Revisions [Date, people, reason]:
%-

base_weight = refiner_weights(2);
d = ones(N+1, 1) * base_weight;
total_tight = 0;
tight_dist = 0;

% The smallest refinement factor allowed
most_refined_weight = refiner_weights(1);

for i=1:size(refinements, 1)
    set = refinements(i, :);
    type = set(3);
    if type == 1
        % In mid cases, we tighten then release after the midpoint
         length = (set(2) - set(1))/2;
         if (mod(length, 1) ~= 0) % Even numbered range pls
             set(2) = set(2) + 1;
             length = (set(2) - set(1))/2;
         end
    else
        length = set(2) - set(1);
    end
    
    if type == 2 && set(2) == N
        set(2) = N+1;
    end
    
    % dx = Br^n
    % Solved for r that makes it reach the most refined at n=length
    r = exp(log(most_refined_weight/base_weight)/length);
    midpoint = false;
    n = 1;
    
    for x=set(1):set(2)
        % Symmetrise the progression for past-mid or min-axial refinements.
        if midpoint == true || type == 0
            d(x) = base_weight + most_refined_weight - base_weight * r^(n-1);
        else
            d(x) = base_weight * r^n;
        end
        
        if d(x) > base_weight
            d(x) = base_weight;
        elseif d(x) < most_refined_weight
            d(x) = most_refined_weight;
        end
        
        n = n + 1;
        if n > length && type == 1 && midpoint == false
            midpoint = true; 
            n = 1;
        end
    end
end
    
summate = sum(d);
factor = L/summate;
d = d * factor;

end

