function x = forward_substitution(A,b, divide)
if nargin == 2
    divide = true;
end

[~,n] = size(A);

x = b;

for i = 1:n
    for j = 1:i-1
        x(i) = x(i) - A(i,j)*x(j);
        add_flops(1); % Just the one multiplication
    end
    if divide
        x(i) = x(i) / A(i,i);
        add_flops(flops_div);
    end
end