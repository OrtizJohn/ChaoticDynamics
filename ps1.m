
%Function is used to control m(number of iterates), x_0 , and R
%respectively
function ps1()
    runMiterates(50,.8,2.5)
end

%Function for logistic map will use variables for each x_n value will be
%represented through the alphabet a,b,c,... respectively with x_0, x_1, and
%so on, R will still represent a constant value that can be changed
function b= logisticalMap(a,R)
    %fprintf('a - %f R - %f ',a,R)
    b = R * a*(1-a);
end
function b= logisticalMap_composition(a,R)
    %fprintf('a - %f R - %f ',a,R)
    b = logisticalMap(logisticalMap(a,R),R);
end
function scatterA(nList,x_nList)
    scatter(nList,x_nList);
    title("x_n versus n ");
    ylim([0 1]);
    xlabel("n");
    ylabel("x_n");
end

function scatterB(x_nList, x_n1List)
    scatter(x_nList,x_n1List);
    title("x_n versus x_n+1 ");
    ylim([0 1]);
    xlim([0 1]);
    xlabel("x_n");
    ylabel("x_n+1");
end
function scatterC(x_nList, x_n2List)
    scatter(x_nList,x_n2List);
    title("x_n versus x_n+2 ");
    ylim([0 1]);
    xlim([0 1]);
    xlabel("x_n");
    ylabel("x_n+2");
end
function runMiterates(m,a,R)
    nList = zeros(1,m);
    x_nList= zeros(1,m);
    x_n1List= zeros(1,m);
    x_n2List= zeros(1,m);
    x_nList(1)= a;
    
    for i=1:m-1
        nList(i+1)= i;
        a= logisticalMap(a,R);
        b=a;
        b=logisticalMap_composition(b,R);
        x_nList(i+1)= a;
        x_n1List(i)= a;
        %fprintf('x_%d -- %f \n',i,a)
        x_n2List(i)= b;
    end
    x_n1List(m) =logisticalMap(a,R);
    x_n2List(m)=logisticalMap_composition(b,R);
    fprintf('\n  n(%d): ',size(nList,2));
    fprintf('%d ', nList);
    fprintf('\n  x_n(%d): ',size(x_nList,2));
    fprintf('%f ', x_nList);
    fprintf('\n  x_n+1(%d): ',size(x_n1List,2));
    fprintf('%f ', x_n1List);
    fprintf('\n  x_n+2(%d): ',size(x_n2List,2));
    fprintf('%f ', x_n2List);
    scatterA(nList,x_nList);
    figure;
    scatterB(x_nList,x_n1List);
    figure;
    scatterC(x_nList,x_n2List);
    
end