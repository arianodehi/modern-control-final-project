clc;
clear;
close all;

R = 6371;
k = 1.2;
r = 0.5;
m = 800;
f = 10000;
g = 9.8;
% Define the differential equations
%space_craft = @(t,x)[
%    dx(1)=x(2)
%    dx(2)=f/m -  g * (R/(R + x(1))) - (k/m)*(x(2))^2 * exp(-x(1)/r);
%    ];
% Initial conditions
X0 = [0 , 0];
time_span = [0,100];
% Solve the differential equations numerically
%[tValues, xValues] = ode45(space_craft, time_span, X0);
% Extract the individual variables
%yValues = xValues(:, 1);
%vValues = xValues(:, 2);
% Plot the output
%figure;
%plot(tValues, yValues, 'r', 'LineWidth', 2);
%hold on;
%plot(tValues, vValues, 'g', 'LineWidth', 2);
%legend('hight' , 'velosity');
%xlabel('Time (seconds)');
%ylabel('Value');

A = [0 1; g*(R/(R+X0(1)))+(k/m)*(X0(2))^2*exp(-X0(1)/r) -2*(k/m)*(X0(2))*exp(-X0(1)/r)];
B = [1;0];
C = [1 0];
D = 0;
X0 = [0 , 0];
%[u , time_span] = gensig('square' , 1000);
%u = heaviside(time_span) * 1000;
%sapac_craft_1 = ss(A,B,C,D);
%[y,t,x] = lsim(sapac_craft_1,u , time_span);
%plot(t,x(:,1),'k',t,x(:,2),'k-.');
%legend('hight' , 'velosity');
%xlabel('Time (seconds)');
%ylabel('Value');
syms 's';
sys = C * (s*eye(2)-A)' * B
%disp(sys);
X0 = [0 , 0];
A = [0 1; g*(R/(R+X0(1)))+(k/m)*(X0(2))^2*exp(-X0(1)/r) -2*(k/m)*(X0(2))*exp(-X0(1)/r)];
B = [1;0];
C = [1 0];
D = 0;
syms 's';
phi  = ilaplace(C * (s * eye(2) - A) * B);
O = ctrb(A,B);
rank(O);
eig(A);
sys = ss2tf(A,B,C,D);
pd = [-3-3i -3+3i];
k1 = place(A,B,pd);
Q = diag([2,5]);
R = 1;
k2 = lqr(A,B,Q,R);
Acl = A-B*k2;
new_sys = ss(Acl,B,C,D);
step(new_sys)
R = 6371;
k = 1.2;
r = 0.5;
m = 800;
f = 10000;
g = 9.8;
X0 = [0 , 0];
A = [0 1; g*(R/(R+X0(1)))+(k/m)*(X0(2))^2*exp(-X0(1)/r) -2*(k/m)*(X0(2))*exp(-X0(1)/r)];
B = [1;0];
C = [1 0];
D = 0;
syms s;
jordan(A);
ctrb(A,B);
obsv(A,C);
eig(A)
x = (3.1305 * eye(2) - A);
rref(x); 
A1 = transpose(A);
C1 = transpose(C);
pd2 = [15-15*i 15+15*i];
g = place(A1,C1,pd2);

