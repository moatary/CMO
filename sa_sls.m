
function  [K3 , vK3, time1, solution, minimum ]=sa_sls(func,itr,subitr)
addpath 'my-optimizers\supp_functs'
evalc(['fun = @', func]);
global all_data
global all_cost
global max_data
global num_of_individuals
global all_scores
global scores
global sel1
global sel2
global xmin
global xmax
[~, xmin, xmax, solution, minimum] = feval(fun);
Nparam=numel(xmin);
xmin(xmin==-Inf)=-10000;
xmax(xmax==Inf)=10000;
X_Max_PSS=xmax; %???
%num_of_individuals=subitr;
max_iteration=itr;
format long

%% definitions 4 SLS
THRESH4SCORE=0.9  ; %%% (TODO) HOW TO SET?
THRESH0=1.5; %%% (TODO)
num_of_individuals=subitr/2; 
max_data=ceil(num_of_individuals*THRESH0);
all_data=zeros(num_of_individuals,Nparam);
all_cost=zeros(num_of_individuals,1);
scores=ones(num_of_individuals,1);
all_scores=ones(max_data,num_of_individuals);
X_Min_PSS=xmin; %???

%% initial
for i=1:num_of_individuals
    rnd=rand(1,numel(xmin(1,:)));
    solution(i,1:Nparam)=  rnd.*(xmax(1,:)-xmin(1,:))+xmin(1,:);
    cost(i,1)=fun( solution(i,1:Nparam));
    if isnan(cost(i,1)) , cost(i,1)=1000000; end
    update_alldata( solution(i,1:Nparam),cost(i,1) ,i)

end

ipop = solution;
ipop_cost = [solution,cost];
ipop_cost_sort = sortrows(ipop_cost,Nparam+1);
G_best = ipop_cost_sort(1,1:Nparam);
G_best_value = ipop_cost_sort(1,Nparam+1);
Tr = [G_best,G_best_value];
  
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
        T = alpha * T;
        
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
        if isnan(y1) , y1=1000000; end
        update_alldata( x1,y1 ,counter);
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
     %% job 4 sls
    for i = 1:num_of_individuals

        LS_new2=update_linear();
 

        costnew2 = fun(LS_new2);
        if isnan(costnew2) , costnew2=1000000; end

        if costnew2 < cost(sel1,1) 
            scores(sel1)=scores(sel1)/THRESH4SCORE;
            all_scores(:,sel1)=1;
          if sel2~=0
              all_scores(sel2,sel1)=all_scores(sel2,sel1)/THRESH4SCORE; 
          end
            
            solution(sel1,1:Nparam) =LS_new2;
            cost(sel1,1) = costnew2;
            update_alldata(LS_new2,costnew2 ,sel1)
           
        else
            update_alldata( LS_new2, costnew2 ,0)
            scores(sel1)=scores(sel1)*THRESH4SCORE;
             if sel2~=0
            all_scores(sel2,sel1)=all_scores(sel2,sel1)*THRESH4SCORE;
             end
       end
    end

    
    ipop_cost = [solution,cost];
    ipop_cost_sort=sortrows(ipop_cost,Nparam+1);
    if ipop_cost_sort(1,Nparam+1) < G_best_value
        G_best = ipop_cost_sort(1,1:Nparam);
        G_best_value = ipop_cost_sort(1,Nparam+1);
        % set the best solution to the one in SA
        x0= G_best;
        y0=G_best_value;
    end
    
    clc;
    K3{Iter} = G_best;
    vK3(Iter) = G_best_value;
    time1 (Iter)=toc;
end
 [~, xmin, xmax, solution, minimum] = feval(fun);
