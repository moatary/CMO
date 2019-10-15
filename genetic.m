function  [K3 , vK3, time1, solution, minimum ]=genetic(func,itr,subitr)
 tic
addpath 'my-optimizers\supp_functs'
 evalc(['constfunc = @', func]);
 [~, xmin, xmax, solution, minimum] = feval(constfunc);
 xmin(xmin==-Inf)=-1000;xmax(xmax==Inf)=1000;
% lb=0.001*ones(length(chosen_subs),1);
% Aeq=[ones(1,sum(chosen_indx==3)),zeros(1,sum(chosen_indx==1));zeros(1,sum(chosen_indx==3)),ones(1,sum(chosen_indx==1))];
% Beq=[1;1];
    options = gaoptimset('Generations',itr,'PopulationSize',subitr);
 vy = ga(@(x)constfunc(x),length(xmin),[],[],[],[],xmin,xmax,[],options);
K3{itr}=vy;
vK3(itr)=constfunc(vy);
time1(itr)=toc;

end
