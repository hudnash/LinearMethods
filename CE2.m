# CE 2 on Linear Methods
## 1. Write a function for multiplying two arbitrary-sized matrices.
### Code Snippet:
```matlab
% Problem #1: Elementary function equivalent to "a*b"
[3 1;8 6;0 4]*[5 9;7 2]
matmult([3 1;8 6;0 4],[5 9;7 2])

function prod = matmult(a,b)
% # columns of the first matrix must equal # rows of the second, where the
% first is a_ik and the second is b_kj. c_ij = sum_k^Nlarger{a_ik*b_kj}.
[ia ja] = size(a); [ib jb] = size(b);
if ja ~= ib
    error("# columns of the first matrix must equal # rows of the second")
end
if ia > jb
    smaller = ib;
else % because we're forcing ja to equal ib, elseif is not necessary. Besides, makes no difference if jb == ia.
    smaller = ja;
end
prod = zeros(ia,jb);
for i=1:1:ia
    for j=1:1:jb
        for k=1:1:smaller
            prod(i,j) = prod(i,j) + a(i,k)*b(k,j);
        end
    end
end
end
```
### Output:
```
ans =

    22    29
    82    84
    28     8


ans =

    22    29
    82    84
    28     8
```
## 2. Solve for CO concentration in all rooms.
### Code Snippet:
```matlab
% Problem #2: Carbon monoxide in a barbecue restaurant controlling levels
% in a children's play area.

% C-O Mass Balances (by room):
% Room 1: 2400+25c3-225c1=0(s.s.)
% Room 2: 225c4-175c2+100=0(s.s.)
% Room 3: 2200+25c1+50c4-275c3=0(s.s.)
% Room 4: 250c3-375c4+25c2=0(s.s.)

% Thus, C.O. concentration in each room computed by 
% [-2400;-100;-2200;0]=[-225 0 25 0;0 -175 0 225;25 0 -275 50;0 25 250
% -375]*[c1;c2;c3;c4]

% Then, conc. = [-225...]\[2400...]:
c = [-225 0 25 0;0 -175 0 225;25 0 -275 50;0 25 250 -375]\[-2400;-100;-2200;0]
```
### Output:
```
c =

   11.8301
   10.4418
   10.4713
    7.6770
```
## 3. Use the inverse matrix to figure out how to reduce CO concentration to less than 10 mg/m^3 in Room 2 (kids section).
### Code Snippet:
```matlab
% And to set c2 = 10 mg/m^3, adjust the inverse matrix [-2400...]^(-1)
% such that:
rat = 1.00;
while c(2) > 10
    c = [-225 0 25 0;0 -175 0 225;25 0 -275 50;0 25 250 -375]\([-2400;-100;-2200;0].*rat);
    rat = rat - .01;
end
disp("Change smoker, room #2, and grill loads to match (respectively):");
[-2400;-100;-2200;0].*rat
```
### Output:
```
Change smoker, room #2, and grill loads to match (respectively):

ans =

       -2256
         -94
       -2068
           0
```
## 4. Write naivegauss function for Gaussian elimination.
### Code Snippet:
```matlab
% Problem #4: Gaussian elimination (naive) method:
function x = naivegauss(a,c)
% Create extended matrix:
a(:,size(a,2)+1)=c; % where c is a column vector!!!
% Forward elimination-- a two-part algorithm loop:
ni = size(a,1); nj = size(a,2);
for i=1:1:ni
% Part I: R1 = R1/a_11
    a(i,:) = a(i,:)/a(i,i);
    % Part II: (R2 --> Rn) = (R2 --> Rn)-(a21 --> an1)*R1
    for k=i+1:1:ni
        a(k,:)=a(k,:)-a(k,i)*a(i,:);
    end
end
% Backwards propagation-- compute each x-value bottom -> top:

x = zeros(size(a,1),1);
% Start the propagation:
lx = length(x);
x(end) = a(lx,lx+1)./a(lx,lx);
for i=lx-1:-1:1
    x(i) = a(i,size(a,2));
    for k=lx:-1:i+1
        x(i) = x(i)-x(k)*a(i,k);
    end
end
end
```
## 5. What are the advantages of LU decomposition over Gaussian elimination?
While Gaussian elimination can run into problems when terms get close to causing NaN and Inf products while forming the upper triangular form, 
LU decomposition avoids the error by forming the resolved matrix from a matrix product. Two matrices, L and U, are used and these reduce the memory
and time spent on larger matrices-- especially when the Gaussian algorithm would otherwise run into errors. LU decomposition utilizes these two factor matrices,
which may be used to reduce unuseful memory and processing time in later matrix operations.

## 6. Find the flow rates of benzene and toluene coming out of the column.
### Code Snippet:
```matlab
% Solved in console, TODO: write new script!!!
```

## 7. Bigger column problem.
### Code Snippet:
```matlab
% Solved in console, TODO: write new script!!!
```

## 8. Write a function to solve a tridiagonal linear system.
### Code Snippet:
```matlab
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
```

## 9. Write a function to multiply two tridiagonal matrices.
### Code Snippet:
```matlab
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
```
