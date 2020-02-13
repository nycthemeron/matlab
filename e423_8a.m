%{
some context: This was for an assignment to find natural frequencies and
mode shapes of a 2DOF system.
%}

clear; close all;
Ms=21;
L=12;
T=1500;

n=4;

m=Ms/n;
l=L/n;

M=zeros(n);
for i=1:n
    M(i,i)=m;
end
K=zeros(n);
for i=1:n
    for j=1:n
        if i==1 && j==1 || i==n && j==n
            K(i,j)=3*T/l;
        elseif i==j
            K(i,j)=2*T/l;
        elseif abs(i-j)==1
            K(i,j)=-T/l;
        end
    end
end

[V D]=eig(K,M)

%Normalized eigenvectors
VNormalized=[V(:,1)*3.5074, V(:,2)*4.5826, V(:,3)*3.5074, V(:,4)*4.5826]

%Adding zeros for endpoints
VN=[0,0,0,0;VNormalized;0,0,0,0];

%This is where the masses are positioned
points=[0, l/2, 3*l/2, 5*l/2, 7*l/2, 12];

%This is for the continuous model
x=0:.1:12;

figure(1)
hold on

%Lumped mass
%The constants VN is multiplied by whatever number is needed to set the
%greatest value to 1. These were calculated before
plot(points,VN(:,1),'-o')
plot(points,VN(:,2),'-o')
plot(points,VN(:,3),'-o')
plot(points,VN(:,4),'-o')

%Continuous
plot(x,sin(pi*x/12))
plot(x,sin(-pi*x/6))
plot(x,sin(-pi*x/4))
plot(x,sin(-pi*x/3))

legend('1st mode','2nd mode','3rd mode','4th mode')
xlabel('position of mass on string')
ylabel('vertical displacement')
title('c. Displacement Profiles for n=4')

%That was messy so I'll divide it into separate plots.
figure(2)
hold on
plot(points,VN(:,1),'-o')
plot(points,VN(:,2),'-o')
plot(points,VN(:,3),'-o')
plot(points,VN(:,4),'-o')

legend('1st mode','2nd mode','3rd mode','4th mode')
xlabel('position of mass on string')
ylabel('vertical displacement')
title('c. Displacement Profiles for n=4, lumped only')

figure(3)
hold on
plot(points,VN(:,1),'-o')
plot(x,sin(pi*x/12))
E1=0;
for i=1:6
    E1=E1+(L/n)*(VN(i,1)-sin(pi*points(i)/12))^2;
end
E1

figure(4)
hold on
plot(points,VN(:,2),'-o')
plot(x,sin(-pi*x/6))
E2=0;
for i=1:6
    E2=E2+(L/n)*(VN(i,2)-sin(pi*points(i)/6))^2;
end
E2

figure(5)
hold on
plot(points,VN(:,3),'-o')
plot(x,sin(-pi*x/4))
E3=0;
for i=1:6
    E3=E3+(L/n)*(VN(i,3)-sin(pi*points(i)/4))^2;
end
E3

figure(6)
hold on
plot(points,VN(:,4),'-o')
plot(x,sin(-pi*x/3))
E4=0;
for i=1:6
    E4=E4+(L/n)*(VN(i,4)-sin(pi*points(i)/3))^2;
end
E4
