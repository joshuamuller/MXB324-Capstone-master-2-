function x = backward_substitution_packed(A,b)
[n,~] = size(A);

x = b;

p = n;
[n,~] = size(b);

for j = n:-1:1
    x(j) = x(j) / A(p);
    p = p - 1;
    add_flops(flops_div);
    for i = j-1:-1:1
        x(i) = x(i) - A(p) * x(j);
        p = p - 1;
        add_flops(1);
    end
end
