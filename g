clc
clear all
clear figure
c = [5 4];   % Objective coefficients

A = [6 4;
    1 2;
    -1 1;
    0 1];  % Converted constraints

b = [24;
    6;
    1;
    2];

% Objective function
Z = @(x1,x2) 5*x1 + 4*x2;

% Constraint functions
c1 = @(x1,x2) 6*x1 + 4*x2 - 24;
c2 = @(x1,x2) 1*x1 +2*x2 - 6;
c3 = @(x1,x2) -1*x1 + 1*x2 -1;
c4 = @(x1,x2) 0*x1 + 1*x2 -2;


%n=size(A, 2);
[m,n] = size(A);
%Phase2 Plotting
%x1=0:max(b./A(:,1)) % inclusive both sides
x1_max = max(b(A(:,1) > 0) ./ A(A(:,1) > 0, 1));
x1 = 0:0.1:x1_max;

for i = 1 : m
    x2 = (b(i) - A(i,1)*x1)/A(i,2)
    plot(x1,x2)
    hold on;
end
%Phase3 All intersection and Corner Points
A = [A; eye(2)];
b = [b; zeros(2,1)];
[m,n] = size(A);
pt = [];
for i = 1 : m
    for j = i + 1 : m
        aa = [A(i,:);A(j,:)];
        bb = [b(i);b(j)];
        if(det(aa) ~= 0)
            X = inv(aa)*bb %aa\bb
            if(X >= 0)
                pt = [pt X]
            end
        end
    end
end
pt=unique(pt','rows')'
%Phase4 Find Feasible Points
FP=[];
z=[];
for i = 1 : size(pt, 2)
    PT1 = pt(1, i);
    PT2 = pt(2, i);
    if(c1(PT1,PT2) <= 0 && c2(PT1,PT2) <= 0 && c3(PT1,PT2) <= 0 && c4(PT1,PT2) <= 0)
        FP = [FP pt(:, i)];
        plot(PT1,PT2,'*r', 'MarkerSize', 10);
        cost = Z(PT1,PT2);
        z = [z cost];
    end
end
disp('FP = ')
disp(FP)
disp('Z = ')
disp(Z)
hold off
%Phase5 Optimal Solution and Optimal Value
[optimal_value,index]=max(z)
optimal_sol=FP(:,index)
