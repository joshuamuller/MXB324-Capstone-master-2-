function x = backward_substitution(A,b, divide)

if nargin == 2
    divide = true;
end

[~,n] = size(A);

x = b;

for i = n:-1:1    % Difference here
    for j = i+1:n % and here
        x(i) = x(i) - A(i,j)*x(j);
        add_flops(1) % A multiplication.
    end
    if divide
        x(i) = x(i) / A(i,i);
        add_flops(flops_div);
    end
end