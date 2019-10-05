function [b, ii] = csr_mult(A, IA, JA, x, rows, columns)
%[b, ii] = csr_mult(A, IA, JA, x) performs matrix by vector multiplication
%          on a matrix in compressed sparse row (CSR) storage. 
%          b = the resulting vector.
%          ii = the diagonal of the original matrix for use elsewhere.

if ~exist('rows', 'var')
    rows = 1:size(x, 1);
end
if ~exist('columns', 'var')
    columns = 1:size(x, 1); 
end
if size(rows', 1) == 0 || size(columns', 1) == 0
    b = 0;
    ii = 1;
    return;
end

b = zeros(size(rows', 1), 1);
ii = b;

index = 1;

for row = 1:rows(1)-1
    amount_in_row = IA(row + 1) - IA(row);
    index = index + amount_in_row;
end

for row=rows
    amount_in_row = IA(row + 1) - IA(row);
    range = index:index+amount_in_row-1;
    i = find(JA(range) == row);
    b(row - rows(1) + 1) = A(range)' * (x(JA(range)) .* ismember(JA(range), columns));
    add_flops(size(JA(range), 1));
    add_flops(size(A(range), 1));
    ii(row - rows(1) + 1) = A(index + i - 1);
    index = index + amount_in_row;
end

end

