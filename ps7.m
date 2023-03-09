function ps7()
    %initial conditions
    %x =  [0 1 2 1 0 0 0 1 0 0 0 1]; % a)
    %x =   [10 -5 2 1 0 0 0 1 0 0 0 1]; % b)
    x =  [0 -1 2 1 0 0 0 1 0 0 0 1]; % c)
    f= lorenz_variational();
    variational_eqns(0,.001,100,x,f);
end


function f = lorenz_variational()
    a = 16;
    r = 45;
    b = 4;
    f=@(x) [a*(x(2)-x(1)), ...              %f1
        (r*x(1))-x(2)-(x(1)*x(3)),...       %f2
        (x(1)*x(2))-(b*x(3)),...            %f3
        a*(x(5)-x(4)),...                   %f4
        (r-x(3))*x(4)-x(5)-x(1)*x(6)...     %f5
        x(2)*x(4)+x(1)*x(5)-b*x(6),...      %f6
        a*(x(8)-x(7)),...                   %f7
        (r-x(3))*x(7)-x(8)-x(1)*x(9),...    %f8
        x(2)*x(7)+x(1)*x(8)-b*x(9),...      %f9
        a*(x(11)-x(10)),...                 %f10
        (r-x(3))*x(10)-x(11)-x(1)*x(12),... %f11
        x(2)*x(10)+x(1)*x(11)-b*x(12)]; %f12

    %f(1) = a*(x(2) - x(1));
    %f(2) = r*x(1) - x(2) - x(1)*x(3);
    %f(3) = x(1)*x(2) - b*x(3);

    %f(4) = a*(x(5)-x(4));
    %f(5) = (r-x(3))*x(4)-x(5)-x(1)*x(6);
    %f(6) = x(2)*x(4)+x(1)*x(5)-b*x(6);
    %f(7) = a*(x(8)-x(7));
    %f(8) = (r-x(3))*x(7)-x(8)-x(1)*x(9);
    %f(9) = x(2)*x(7)+x(1)*x(8)-b*x(9);
    %f(10) = a*(x(11)-x(10));
    %f(11) = (r-x(3))*x(10)-x(11)-x(1)*x(12);
    %f(12) = x(2)*x(10)+x(1)*x(11)-b*x(12);
    
end 

function variational_eqns(t0, h, n, x, F)
    %x is augmented state vector [x y z &xx &xy &xz &yx &yy &yz &zx %&zy &zz]
    t = t0;
    w= zeros(n+1,length(x));
    w(1,:) = x;
    ct=0;
    for i = 1:n+1
        
        fprintf('%d - State vector at time %f: %.10f , %.10f , %.10f , %.10f , %.10f , %.10f , %.10f , %.10f , %.10f , %.10f , %.10f , %.10f \n \n',ct, t, w(i,1), w(i,2), w(i,3), w(i,4), w(i,5), w(i,6), w(i,7), w(i,8), w(i,9), w(i,10), w(i,11), w(i,12));   
        %RK4
        ct=ct+1;
        x =w(i,:);
        K1 = F(x);
        K2 = F(w(i,:) + h/2*K1);
        K3 = F(w(i,:) + h/2*K2);
        K4 = F(w(i,:) + h*K3);
        
        w(i+1,:) = w(i,:) + h*((K1 + 2*K2 + 2*K3 + K4)/6);
        
        t = t + h; 
    end
    X = w(i+1,:);
    ev1 = X(4) + X(5) + X(6);
    ev2 = X(7) + X(8) + X(9);
    ev3 = X(10) + X(11) + X(12);
    fprintf('Evolved variations: \n %.10f , %.10f , %.10f ', ev1,ev2,ev3);

end
