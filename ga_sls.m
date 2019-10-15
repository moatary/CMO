function [K3 , vK3, time1, solution, minimum ]=ga_sls(func,itr,subitr)
    %%% (VIMP)
PRECI = 22;    % Precisicion of binary representation

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
    X_Max_PSS=xmax;
%%
THRESH4SCORE=0.9  ; %%% (TODO) HOW TO SET?
THRESH0=1.5; %%% (TODO)
num_of_individuals=subitr/2; 
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
    Iop_PSS(i,1:Nparam)=  rnd.*(xmin(1,:)-xmin(1,:))+(1-rnd).*xmin(1,:);
    Cost(i,1)=fun( Iop_PSS(i,1:Nparam));
    update_alldata( Iop_PSS(i,1:Nparam),Cost(i,1) ,i)

end

ipop_cost = [Iop_PSS,Cost];
%%% (TODO) Ino ke baese complexity mishe hazf kon va tu for emalesh kon
ipop_cost_sort = sortrows(ipop_cost,Nparam+1);
G_best = ipop_cost_sort(1,1:Nparam);
G_best_value = ipop_cost_sort(1,Nparam+1);
Tr = [G_best,G_best_value];

%%


for Iter=1:itr
    tic
    
      NIND = num_of_individuals; % Number of individuals per subpopulations
    GGAP = 1;           % Generation gap, how many new individuals are created
    SEL_F = 'rws';       % Name of selection function
    XOV_F = 'xovsp';     % Name of recombination function for individuals
    MUT_F = 'mut';       % Name of mutation function for individuals
    FieldDR =[X_Min_PSS;X_Max_PSS];% feval(OBJ_F,[],1); %%%
    NVAR = size(FieldDR,2);
    %%% first making iop binary
    sgn=sign(Iop_PSS);
    bins=dec2bin(floor(abs(Iop_PSS(:)).*(2^(PRECI-6))));
    PRECI2=size(bins,2);
    bins1=arrayfun(@(n)str2num(transpose(bins(n,:))),1:size(bins,1),'UniformOutput', false);
    bins2=cell2mat(bins1);bins2=transpose(bins2);
    Chrom=zeros(num_of_individuals,Nparam*PRECI2);
    for i=1:Nparam
        Chrom(:,(i-1)*PRECI2+1:i*PRECI2)=bins2((i-1)*num_of_individuals+1:i*num_of_individuals,:);
    end
    ObjV=Cost(:,1);
    FieldDD = [rep([PRECI2],[1, NVAR]);...
        FieldDR;...
        rep([1; 0; 1 ;1], [1, NVAR])];
    FitnV = ranking(ObjV);
    [ SelCh ,~]= select(SEL_F, Chrom, FitnV, GGAP);
    SelCh=recombin(XOV_F, SelCh);
    SelCh=mutate(MUT_F, SelCh);
    matbin=repmat([2.^(PRECI2-1:-1:0)],size(Iop_PSS,1),size(Iop_PSS,2));
    matbin=matbin.*SelCh;
    matbin2=arrayfun(@(n)(sum(matbin(:,(n-1)*PRECI2+1:n*PRECI2),2))./(2^(PRECI-6)),1:size(Iop_PSS,2),'UniformOutput' , false);
    TL_new_PSS=cell2mat(matbin2).*sgn;
    cost_new=zeros(size(   TL_new_PSS,1) ,1);
    for i=1:size(   TL_new_PSS,1)
        cost_new(i,1) = fun(TL_new_PSS(i,1:Nparam));
    end
    costpss2=[cost_new,TL_new_PSS];
    costpss2=sortrows(costpss2,1);
    for m=1:numel(cost_new)
        [cs,ind]=max(Cost);
        if costpss2(m,1)<cs(1)
            Cost(ind(1))= costpss2(m,1);
            Iop_PSS(ind(1),:)=costpss2(m,2:end);
            update_alldata(costpss2(m,2:end),costpss2(m,1) ,ind(1))

        else
            break
        end
    end
    
    ipop_cost = [Iop_PSS,cost_new];
    ipop_cost_sort=sortrows(ipop_cost,Nparam+1);
    if ipop_cost_sort(1,Nparam+1) < G_best_value
        G_best = ipop_cost_sort(1,1:Nparam);
        G_best_value = ipop_cost_sort(1,Nparam+1);
    end

 
    
    
    
    %% job 4 sls
    for i = 1:num_of_individuals

        LS_new2=update_linear();
 

        costnew2 = fun(LS_new2);

        if costnew2 < Cost(sel1,1)
            scores(sel1)=scores(sel1)/THRESH4SCORE;
            all_scores(:,sel1)=1;
          if sel2~=0
              all_scores(sel2,sel1)=all_scores(sel2,sel1)/THRESH4SCORE; 
          end
            
            Iop_PSS(sel1,1:Nparam) =LS_new2;
            Cost(sel1,1) = costnew2;
            update_alldata(LS_new2,costnew2 ,sel1)
           
       else
            update_alldata( LS_new2, costnew2 ,0)
            scores(sel1)=scores(sel1)*THRESH4SCORE;
             if sel2~=0
            all_scores(sel2,sel1)=all_scores(sel2,sel1)*THRESH4SCORE;
             end
       end
    end

    
    ipop_cost = [Iop_PSS,cost_new];
    ipop_cost_sort=sortrows(ipop_cost,Nparam+1);
    if ipop_cost_sort(1,Nparam+1) < G_best_value
        G_best = ipop_cost_sort(1,1:Nparam);
        G_best_value = ipop_cost_sort(1,Nparam+1);
    end
    %Tr = [G_best,G_best_value];
        K3{Iter} = G_best;
vK3(Iter) = G_best_value;
 time1 (Iter)=toc;
end
end
