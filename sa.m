
function  [K3 , vK3, time1, solution, minimum ]=sa(func,itr,subitr)
addpath 'my-optimizers\supp_functs'
evalc(['fun = @', func]);

[Nparam, min_range, max_range, solution, minimum] = feval(fun);
min_range(min_range==-Inf)=-10000; %%%
max_range(min_range==Inf)=10000; %%%
num_of_individuals=subitr;
format long
global xmin
global xmax
[~, xmin, xmax, solution, minimum] = feval(fun);
Nparam=numel(xmin);
xmin(xmin==-Inf)=-10000;
xmax(xmax==Inf)=10000;
%num_of_individuals=subitr;


%% job 4 sa
    rnd=rand(1,numel(xmin(1,:)));
    x0=  rnd.*(xmax(1,:)-xmin(1,:))+xmin(1,:);
    y0 = fun(x0);
x = x0;
y = y0;
alpha = 0.6;
T = 0.1;
s = 0.2* 0.25*(min(abs(xmin))+min(abs(xmax))+max(abs(xmin))+max(abs(xmax)));


accept=0;
j=1;l=1;
for Iter=1:itr
    tic
    
    counter=0;
    
    while 1==1
        counter=counter+1;
        
        if j==21
            j=1;  T = alpha * T;
        end;
        if rem(l,100)==0
            if accept < 25, s = s / 2;
            end
            if accept > 75, s = 2 * s;
            end
            x0 = x;
            y0 = y;
            
            accept=0;
            l=1;
            j=j+1;
        end
        x1 = x0 + (0.5 - rand (size (x0))) * s;
        y1 = fun(x1);
        dy = y1 - y0;
        if rand < exp (- dy / T);
            x0 = x1;
            y0 = y1;
            accept = accept + 1;
            if y0 < y;
                x = x0;
                y = y0;
            end
        end
        l=l+1;
        
        
        if counter>=num_of_individuals,   break, end
    end
  
    clc;
    K3{itr}=x;
    vK3(itr)=y;
    time1(itr)=toc;
end
[~, xmin, xmax, solution, minimum] = feval(fun);
