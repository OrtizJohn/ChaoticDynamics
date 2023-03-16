%4th order runge-kutta 
function ps4()
    clc;
    clear all;
   
    %testLorenz()
    
    %DampPendulum()
    PendulumPortrait()
end

function DampPendulum()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%q1 defining function / run rk4
%   Use m=0.1kg, l=0.1m, β=0 and set the drive amplitude and frequency to zero (α = A = 0).
%   (a) Turn in a plot of the state-space trajectory emanating from the point [θ, ω] = [3, 0.1] with
%   ∆t = 0.005. Is this initial condition near a fixed point? Which one? Is that point stable or
%   unstable?
% θ''(t) = (A/mL)*cos(αt) - (β/m) * θ'(t) - (g/L) sin(θ(t)) 
% note this equation is non-autonomous b/c  there is explicit t on right hand side (so function does depend on time)
% θ'' = (A/mL)*cos(αt) - (β/m) * θ' - (g/L) sin(θ)
% will define θ as θ_1 to not get confused:
% θ_1''(t) = (A/mL)*cos(αt) - (β/m) * θ_1'(t) - (g/L) sin(θ_1(t)) 
% defining n-1 helper variables: θ_1' = θ_2
% Subbing in helper variable: θ_2' = (A/mL)*cos(αt) - (β/m) * θ_2 - (g/L) * sin(θ)
% Getting 2 1st order equations 
% θ_1' = θ_2
% θ_2' = (A/mL)*cos(αt) - (β/m) * θ_2 - (g/L) * sin(θ)
%
%Sanity Check (plug back in )
% θ_1'' = (A/mL)*cos(αt) - (β/m) * θ_1' - (g/L) * sin(θ_1) = θ''(t) = (A/mL)*cos(αt) - (β/m) * θ'(t) - (g/L) sin(θ(t))
% Pass!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defining variables
m =.1; % kg
L = .1; % m 
beta = 0; 
alpha = 0; % drive ampltitude
A = 0; %drive frequency
x = [.01,0]; %initial point
h = .005;
g = 9.8;
f= @(t,x) [x(2),(A*cos(alpha*t) - beta*L*x(2) - m*g*sin(x(1)))/(m*L)];
[t,w] = rk4(0,h,500,x,f);
hold on;
plot(w(:,1),w(:,2),'Marker','.')
title(['State Space for forced, damped pendulum for IC = (' num2str(x(1)) ', ' num2str(x(2)) ')']);
xlabel('theta');
ylabel('w');
end
function PendulumPortrait()
%defining variables
m =.1; % kg
L = .1; % m 
beta = 0; 


%natural frequency = sqrt(g/l)/2pi
g = 9.8;
nf = sqrt(g/L)/(2*pi); %natural frequency 
 
alpha = 0;%.5*nf; % %drive frequency 
A = 0; %ampltitude .75
h = .05;

f= @(t,x) [x(2),(A*cos(alpha*t) - beta*L*x(2) - m*g*sin(x(1)))/(m*L)];

hold on;

a=-10;
b=10;
for i=1:250
    x = [a+(b-a).*rand(1,1), a+(b-a).*rand(1,1)];
    
    [t,w] = rk4(0,h,250,x,f);
    scatter(mod(w(:,1),2*pi),w(:,2),'Marker','.','Color',[.18, .835, .784])
    
end




title(['State Space for forced, damped pendulum Portrait']);
xlabel('theta mod 2pi');
ylabel('w');
%xlim([-10,10])
hold off
end

function testLorenz()
    %variables 
    a=16;
    r=50;
    b=4;
    x= [0,1,2]; %initial condition
    h=.001;
    n=3;


    f=@(t,vec) [a*(vec(2)-vec(1)),vec(1)*r-vec(2)*vec(3)-vec(2),vec(1)*vec(2)-b*vec(3)];

    [t,w] = rk4(0,h,n,x,f);
    hold on;
    wT = transpose(w);
    plot3(wT(1,:),wT(2,:),wT(3,:))
    
end
% starting time t0, time step ∆t, number of steps n, and starting value ~x(t0) for the
%state vector
function [t,w] = rk4(t0,deltaT,n,x,f)
    %fprintf('generating row vector from %f to %f with step size of %f \n',t0,t0+deltaT*n, deltaT);
    t= linspace(t0,t0+deltaT*n,n+1);
    w= zeros(n,length(x));

    w(1,:)= x;
    for i= 1:n
        k1 = deltaT* f(t(i),w(i,:));
        k2 = deltaT*f(t(i)+deltaT/2.0,w(i,:)+k1/2.0);
        k3 = deltaT*f(t(i)+deltaT/2.0,w(i,:)+k2/2.0);
        k4 = deltaT*f(t(i)+deltaT,w(i,:)+k3);
        w(i+1,:) = w(i,:) + (k1 + 2.0*(k2+k3) + k4)/6.0;
        %fprintf('After %d timestep: x:%f y:%f z:%f \n',i,w(i+1,1),w(i+1,2),w(i+1,3)) %for test lorenz
        %fprintf('After %d timestep: theta:%f w:%f  \n',i,w(i+1,1),w(i+1,2)) %for ps4
    end    
end