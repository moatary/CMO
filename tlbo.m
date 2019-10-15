function [K3 , vK3 , time1 , solution, minimum]=tlbo(func,itr,subitr)
addpath 'my-optimizers\supp_functs'
evalc(['fun = @',func]);
[~, xmin, xmax, solution, minimum] = feval(fun);
     xmin(xmin==-Inf)=-1000;xmax(xmax==Inf)=1000;

Nparam=numel(xmin);

X_Min_PSS=xmin;
X_Max_PSS=xmax;
Number_of_Pop=subitr/2;

%% initial
for i=1:Number_of_Pop
    rnd=rand(1,numel(xmin(1,:)));
    Iop_PSS(i,1:Nparam)=  rnd.*(xmin(1,:)-xmin(1,:))+(1-rnd).*xmin(1,:);
    Cost(i,1)=fun( Iop_PSS(i,1:Nparam));
end

ipop = Iop_PSS;
ipop_cost = [Iop_PSS,Cost];
ipop_cost_sort = sortrows(ipop_cost,Nparam+1);
G_best = ipop_cost_sort(1,1:Nparam);
G_best_value = ipop_cost_sort(1,Nparam+1);
Tr = [G_best,G_best_value];


%%
for Iter=1:itr
    tic
    mean1=sum(ipop_cost)./Number_of_Pop;  % (???)
    
    
    
    for i=1:Number_of_Pop
        
        [TL_new_PSS(i,:)]=Teacher_PSS(X_Max_PSS,X_Min_PSS,Nparam,0,mean1,Tr,Iop_PSS(i,:));
        
        Cost_learned(i,1) = fun(TL_new_PSS(i,:));
        if Cost_learned(i,1) < Cost(i,1);
            ipop_new(i,:) = TL_new_PSS(i,:);
            cost_new(i,1) = Cost_learned(i,1);
        else
            ipop_new(i,:) = ipop(i,:);
            cost_new(i,1) = ipop_cost(i,Nparam+1);
        end
    end
    
    for i = 1:Number_of_Pop/2
        LL = randperm(Number_of_Pop);
        if LL(1)==i
            LL(1)=[]; %(???) chera?
        elseif LL(2)==i
            LL(2)=[]; %(???) chera?
        end
        
        [LS_new_PSS(i,:)]=Student_PSS(LL,Iop_PSS(i,:),Iop_PSS,X_Max_PSS,X_Min_PSS,Nparam);
        Cost_learned1(i,1) = fun(LS_new_PSS(i,:));
        
        ipop_new1(i,:) = LS_new_PSS(i,:);
        cost_new1(i,1) = Cost_learned1(i,1);
        if cost_new1(i,1) < cost_new(i,1)
            ipop_new(i,1:Nparam) = ipop_new1(i,1:Nparam);
            cost_new(i,1) = cost_new1(i,1);
        end
    end
    Cost = cost_new;
    Iop_PSS = ipop_new;
    ipop_cost = [ipop_new,cost_new];
    ipop_cost_sort=sortrows(ipop_cost,Nparam+1);
    if ipop_cost_sort(1,Nparam+1) < G_best_value
        G_best = ipop_cost_sort(1,1:Nparam);
        G_best_value = ipop_cost_sort(1,Nparam+1);
    end
    K3{Iter} = G_best;
    vK3(Iter) = G_best_value;
    time1(Iter)=toc;
end

end