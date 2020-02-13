%{
some context: This was written for a lab involving a pendulum on a rotary
beam. In this lab, we found parameters of the system and then used output
feedback control to stabilize the inverted pendulum.
%}

clear
close all
load('C:\Users\redli\Documents\ENAE464 - Aero Lab II\all_data.mat')

%%Taking the omega and sigma values from each data set...

omega1=6.7003;
sigma1=0.0455;

omega2=8.3777;
sigma2=0.8843;

omega3=7.9468;
sigma3=0.5279;

%k and c are "K-phi" and "K-phidot"
k=6;
c=0.8;
%r0 is the distance along beam from pivot to disk weights, in meters
r0=(6+5/8)/39.37;
%m0 is the mass of the disk weights, in kg
m0=65.1/1000;
%r is the distance along beam from pivot to pendulum, in meters
r=(6+5/8)/39.37;
%m is the mass of the pendulum, in kg
m=125.4/1000;
%g is the acceleration due to gravity
g=9.81;
%l is the length of the pendulum, in meters
l=(13+1/8)/39.37;

phi=data4.phi;
theta=data4.theta;
phi_dot=data4.phi_dot;
theta_dot=data4.theta_dot;
time=data4.time;

%Solve for motor constants K1, K2, and inertia J.
M = [-c 1 -2*sigma2;
k 0 -omega2^2-sigma2^2;
-c 1 -2*sigma3^2;
k 0 -omega3^2-sigma3^2];
Y = [0 0 2*sigma3^2*m0*r0^2 (omega3^2+sigma3^2)*m0*r0^2]';
X = M\Y;
K1 = X(1)
K2 = X(2)
J = X(3)

%%Create matrices for model
A = [0 0 1 0; 0 0 0 1;
0 m*g*r/J -K2/J 0; ...
0 -(J+m*r^2)*g/J/l r*K2/J/l 0];
B = [0 0 K1/J -r*K1/J/l]';
C = [1 0 0 0; 0 1 0 0];
D = [0 0]';
K = [k 0 -c 0];

% Define linear system
sys = ss(A-B*K,B,C-D*K,D);

% Specify initial condition of the simulation from experimental data
x0= [phi(1),theta(1),phi_dot(1),theta_dot(1)];

% Run simulation
[y,t] = initial(sys,x0,10);

% Plot results. Compare model to experimental data.
subplot(211)
plot(t,y(:,1));hold on
plot(time-time(1),phi)
legend('Model','Data')
ylabel('Phi'),xlabel('Time')
subplot(212)
plot(t,y(:,2));hold on
plot(time-time(1),theta)
ylabel('Theta'),xlabel('Time')