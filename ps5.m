
function ps5()
    clc;
    clear all;
    view(3)
    %Lorenz();
    Rossler();
end
function Lorenz()
    %variables 
    a=16;
    r=45;
    b=4;
    x= [-13,-12,52]; %initial condition
    h=.01;
    n=500;

    f=@(x) [a*(x(2)-x(1)),(r*x(1))-x(2)-(x(1)*x(3)),(x(1)*x(2))-(b*x(3))]; %lorenz
    
    tol=.01;
    [t,w] = runAdaptiveRk4(0,h,n,x,f,tol);
    %[t2,w2] = rk4(0,.001,n,x,f);


    %plotting results 
    hold on;
    title(['Lorenz system with IC=[-13,-12,52] using AdaptiveRk4 ErrorTolerance: ',string(tol)]); %q2a
    %title('Lorenz system with IC=[-13,-12,52] AdaptiveRk4 vs Rk4 n=1000'); %q2b
    xlabel('X');
    ylabel('Y');
    zlabel('Z'); 
    scatter3(w(:,1),w(:,2),w(:,3),5,"MarkerFaceColor",[0.3010 0.7450 0.9330],"MarkerEdgeColor",[0.3010 0.7450 0.9330],"Marker",".")
    %scatter3(w2(:,1),w2(:,2),w2(:,3),5,"MarkerFaceColor",[0.9290 0.6940 0.1250],"MarkerEdgeColor",[0.9290 0.6940 0.1250],"Marker",".")
    %p4 plot h = .05, .5 n=225,500 tol= 15
end
function Rossler()
    %variables 
    a=.398;
    b=2;
    c=4;
    x= [-13,-12,52]; %initial condition
    h=.01;
    n=500;

    %f=@(x) [a*(x(2)-x(1)),(r*x(1))-x(2)-(x(1)*x(3)),(x(1)*x(2))-(b*x(3))]; %lorenz
    f=@(x) [-1*(x(2) +x(3)), x(1) + a*x(2), b + x(3)*(x(1)-c)];
    tol=.025;
    [t,w] = runAdaptiveRk4(0,h,n,x,f,tol);
    


    %plotting results 
    hold on;
    title(['Rossler system with IC=[-13,-12,52] using AdaptiveRk4 ErrorTolerance: ',string(tol)]); %q2a
    xlabel('X');
    ylabel('Y');
    zlabel('Z'); 
    scatter3(w(:,1),w(:,2),w(:,3),5,"MarkerFaceColor",[0.3010 0.7450 0.9330],"MarkerEdgeColor",[0.3010 0.7450 0.9330],"Marker",".")
   
end
function [t,w] = rk4(t0,deltaT,n,x,f)
    tVect = zeros(n,1);
    w= zeros(n,length(x));
    
    w(1,:)= x;
    i=1;
    t=t0;
    finalT =deltaT*n;
    while t<= finalT
        %fprintf('At %f timestep: x,y,z: %f %f %f \n',t,w(i,1),w(i,2),w(i,3))%for  lorenz
        %time step deltaT with rk4
        k1 = f(w(i,:));
        k2 = f(w(i,:)+ deltaT/2.0*k1);
        k3 = f(w(i,:)+ deltaT/2.0*k2);
        k4 = f(w(i,:)+ deltaT*k3);
        w(i+1,:) = w(i,:) + deltaT*((k1 + 2.0*(k2+k3) + k4)/6.0);
        
        %update new time 
        t = t+deltaT;
        tVect(i+1) = t;
        %update i
        i = i+1;

    end
end
function [t,w] = runAdaptiveRk4(t0,deltaT,n,x,f,tol)
%function will run adaptive 4th order runge kutta using an adaptive timestep
    tVect = zeros(n,1);
    w= zeros(n,length(x));
    
    w(1,:)= x;
    i=1;
    t=t0;
    timestep0= deltaT;
    finalT =timestep0*n;
    while t<= finalT
        fprintf('At %f timestep: x,y,z: %f %f %f -- %f\n',t,w(i,1),w(i,2),w(i,3),deltaT)%for  lorenz
        %time step deltaT with rk4
        k1 = f(w(i,:));
        k2 = f(w(i,:)+ deltaT/2.0*k1);
        k3 = f(w(i,:)+ deltaT/2.0*k2);
        k4 = f(w(i,:)+ deltaT*k3);
        w(i+1,:) = w(i,:) + deltaT*((k1 + 2.0*(k2+k3) + k4)/6.0);

        %time step deltaT/2 
        k1 = f(w(i,:));
        k2 = f(w(i,:)+ deltaT/4.0*k1);
        k3 = f(w(i,:)+ deltaT/4.0*k2);
        k4 = f(w(i,:)+ deltaT/2.0*k3);
        
        tempW  = w(i,:) + deltaT/2.0*((k1 + 2.0*(k2+k3) + k4)/6.0);
        
        k1 = f(w(i,:));
        k2 = f(w(i,:)+ deltaT/4.0*k1);
        k3 = f(w(i,:)+ deltaT/4.0*k2);
        k4 = f(w(i,:)+ deltaT/2.0*k3);
        
        x2  = tempW + deltaT/2.0*((k1 + 2.0*(k2+k3) + k4)/6.0);

        %error
        deltaW = x2 - w(i+1,:) ;
        
        error = max(abs(deltaW));
        
        %update new time 
        t = t+deltaT;
        tVect(i+1) = t;
        %update state vector
        w(i+1,:) = x2;

        

        %update timestep deltaT dependent on error
        if(error == 0)
            deltaT = deltaT;
        elseif( error >= tol)
            deltaT = deltaT*(abs(tol/error))^.2;
        end
        %update i
        i = i+1;

    end    

end
