function x = forward_substitution_packed(A,b)
[n,~] = size(b);

x = b;
p = 1;

for j = 1:n
    for i = 1:j-1
        x(j) = x(j) - A(p)*x(i);
        p = p + 1;
        add_flops(1);
    end
    x(j) = x(j) / A(p);
    p = p + 1;
    add_flops(flops_div);
end