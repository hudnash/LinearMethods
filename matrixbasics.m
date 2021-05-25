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

% And to set c2 = 10 mg/m^3, adjust the inverse matrix [-2400...]^(-1)
% such that:
rat = 1.00;
while c(2) > 10
    c = [-225 0 25 0;0 -175 0 225;25 0 -275 50;0 25 250 -375]\([-2400;-100;-2200;0].*rat);
    rat = rat - .01;
end
disp("Change smoker, room #2, and grill loads to match (respectively):");
[-2400;-100;-2200;0].*rat

% Problem 4: (see end of script) "Naive" Gauss-Jordan elimination:
[-225 0 25 0;0 -175 0 225;25 0 -275 50;0 25 250 -375]\[-2400;-100;-2200;0]
naivegauss([-225 0 25 0;0 -175 0 225;25 0 -275 50;0 25 250 -375],[-2400;-100;-2200;0])

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