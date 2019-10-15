function [K3 , vK3 , time1, solution, minimum ]=sls(func,itr,subitr)
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
     xmin(xmin==-Inf)=-1000;xmax(xmax==Inf)=1000;

Nparam=numel(xmin);
%%
THRESH4SCORE=0.9  ; %%% (TODO) HOW TO SET?
THRESH0=1.5; %%% (TODO)
num_of_individuals=subitr; 
max_data=ceil(num_of_individuals*THRESH0);
all_data=zeros(num_of_individuals,Nparam);
all_cost=zeros(num_of_individuals,1);
scores=ones(num_of_individuals,1);
all_scores=ones(max_data,num_of_individuals);


format long
%itr=15;

X_Min_PSS=xmin;


%% initial
for i=1:num_of_individuals
    rnd=rand(1,numel(xmin(1,:)));
    ipop_new(i,1:Nparam)=  rnd.*(xmin(1,:)-xmin(1,:))+(1-rnd).*xmin(1,:);
    cost_new(i,1)=fun( ipop_new(i,1:Nparam));
    update_alldata( ipop_new(i,1:Nparam),cost_new(i,1) ,i)
end

%ipop = ipop_new;
ipop_cost = [ipop_new,cost_new];
ipop_cost_sort = sortrows(ipop_cost,Nparam+1);
G_best = ipop_cost_sort(1,1:Nparam);
G_best_value = ipop_cost_sort(1,Nparam+1);
%Tr = [G_best,G_best_value];

%%
for Iter=1:itr
    %%
    tic
    for i = 1:num_of_individuals

        LS_new2=update_linear();
 

        costnew2 = fun(LS_new2);

        if costnew2 < cost_new(sel1,1)
            scores(sel1)=scores(sel1)/THRESH4SCORE;
            all_scores(:,sel1)=1;
          if sel2~=0
              all_scores(sel2,sel1)=all_scores(sel2,sel1)/THRESH4SCORE; 
          end
            
            ipop_new(sel1,1:Nparam) =LS_new2;
            cost_new(sel1,1) = costnew2;
            update_alldata(LS_new2,costnew2 ,i)
           
       else
            update_alldata(LS_new2, costnew2 ,0)
            scores(sel1)=scores(sel1)*THRESH4SCORE;
             if sel2~=0
            all_scores(sel2,sel1)=all_scores(sel2,sel1)*THRESH4SCORE;
             end
       end
    end

    
 
    ipop_cost = [ipop_new,cost_new];
    ipop_cost_sort=sortrows(ipop_cost,Nparam+1);
    if ipop_cost_sort(1,Nparam+1) < G_best_value
        G_best = ipop_cost_sort(1,1:Nparam);
        G_best_value = ipop_cost_sort(1,Nparam+1);
    end
  %  Tr = [G_best,G_best_value];
        K3{Iter} = G_best;
vK3(Iter) = G_best_value;
time1(Iter)=toc;
end

end