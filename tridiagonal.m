
% Problem 8: tridiagsolve(a,b,c,d)
A = [2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
a = [2 2 2 2]; b = [-1 -1 -1]; c = [-1 -1 -1];
d = [1 2 3 4]';
A \ d
B = [2 -1 0 0; 0 1 -1 0; 0 0 1.6667 -1; 0 0 0 1.75];
a2=[2 1 1.667 1.75];b2=[0 0 0];c2=[-1 -1 -1];
tridiagsolve(b,a,c,d)
x3 = A*B
x2 = tridiagmult(b,a,c,b2,a2,c2)
function x = tridiagsolve(b,a,c,d)
% Given:[ a1 c1 0  0  d1;   
%         b2 a2 c2 0  d2;
%         0  b3 a3 c3 d3;
%         0  0  b4 a4 d4 ]
b = [0 b]; c = [c 0];
n = length(a);
x = zeros(n,1);
% Forward elimination
for i=1:1:n-1
    % c(i) = c(i)/a(i); % Incorrect, because skip /b(i)
    % a(i) = 1;         % Incorrect, because skip /b(i)
    % d(i) = d(i)/a(i); % Incorrect, because skip /b(i)
    a(i+1) = a(i+1)-c(i)*b(i+1)/a(i);
    d(i+1) = d(i+1)-d(i)*b(i+1)/a(i);
    b(i+1) = 0;
end
% Back substitution
x(end) = d(end)/a(end);
for i=n-1:-1:1
    x(i) = (d(i)-c(i)*x(i+1))/a(i);
end
end

function prod = tridiagmult(a1,b1,c1,a2,b2,c2)
% Given:[ b1 c1 0  0  d1;   
%         a2 b2 c2 0  d2;
%         0  a3 b3 c3 d3;
%         0  0  a4 b4 d4 ]
a1 = [0 a1]; c1 = [c1 0]; a2 = [0 a2]; c2 = [c2 0];
if length(b1) ~= length(b2)
    error("matrices must be the same size square matrices")
end
n = length(b1);
for i=1:1:n
    for j=1:1:n
        if i-1 == j
            prod(i,j)=a1(i)*b2(j);
        elseif i == j
           
        elseif
        elseif
        end
    end
end
end