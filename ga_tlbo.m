function  [K3 , vK3, time1, solution, minimum ]=ga_tlbo(func,itr,subitr)
addpath 'my-optimizers\supp_functs'
evalc(['fun = @',func]);
    [~, xmin, xmax, solution, minimum] = feval(fun);
         xmin(xmin==-Inf)=-1000;xmax(xmax==Inf)=1000;

    Nparam=numel(xmin);
X_Min_PSS=xmin;
X_Max_PSS=xmax;
num_of_individuals=subitr/2;

%% initial
for i=1:num_of_individuals
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
    mean1=sum(ipop_cost)./num_of_individuals;  % (???)
    
    for i=1:num_of_individuals
        
     [TL_new_PSS(i,:)]=Teacher_PSS1(X_Max_PSS,X_Min_PSS,Nparam,0,mean1,Tr,Iop_PSS(i,:));
 
Cost_learned(i,1) = fun(TL_new_PSS(i,:));
     
        if Cost_learned(i,1) < Cost(i,1);    %%% (INJA DARE ETELAAT HAZF MIKONE)
            ipop_new(i,:) = TL_new_PSS(i,:);
            cost_new(i,1) = Cost_learned(i,1);
        else
            ipop_new(i,:) = ipop(i,:);
            cost_new(i,1) = ipop_cost(i,Nparam+1);
        end
    end
    
    
    
    
    
    %% STUDENT PHASE:
    
    
NIND = num_of_individuals; % Number of individuals per subpopulations
GGAP = 1;           % Generation gap, how many new individuals are created
SEL_F = 'rws';       % Name of selection function
XOV_F = 'xovsp';     % Name of recombination function for individuals
MUT_F = 'mut';       % Name of mutation function for individuals

% Get boundaries of objective function
   FieldDR =[X_Min_PSS;X_Max_PSS];% feval(OBJ_F,[],1); %%%

% Number of variables of objective function, in OBJ_F defined
   NVAR = size(FieldDR,2);   

% Build fielddescription matrix
   PRECI = 22;    % Precisicion of binary representation
 %  prec10=floor(log(2^PRECI)/log(10));
   

%%% first making iop binary
sgn=sign(ipop_new);
    bins=dec2bin(floor(abs(ipop_new(:)).*(2^(PRECI-6))));
PRECI2=size(bins,2);
bins1=arrayfun(@(n)str2num(transpose(bins(n,:))),1:size(bins,1),'UniformOutput', false);
bins2=cell2mat(bins1);bins2=transpose(bins2);
Chrom=zeros(num_of_individuals,Nparam*PRECI2);
for i=1:Nparam
    Chrom(:,(i-1)*PRECI2+1:i*PRECI2)=bins2((i-1)*num_of_individuals+1:i*num_of_individuals,:);
end
ObjV=cost_new(:,1);


   FieldDD = [rep([PRECI2],[1, NVAR]);...
              FieldDR;...
              rep([1; 0; 1 ;1], [1, NVAR])];
          
          
 % then all crossover stuff
%%% +(TODO) some mods

   % Fitness assignement to whole population
      FitnV = ranking(ObjV);
            
   % Select individuals from population
     [ SelCh ,~]= select(SEL_F, Chrom, FitnV, GGAP);
     
   % Recombine selected individuals (crossover)
      SelCh=recombin(XOV_F, SelCh);

   % Mutate offspring
      SelCh=mutate(MUT_F, SelCh);

   % Insert offspring in population replacing parents
   %   Chrom = reins(Chrom, SelCh);

% then evaluate again each values % reset count variables


matbin=repmat([2.^(PRECI2-1:-1:0)],size(ipop_new,1),size(ipop_new,2));
matbin=matbin.*SelCh;
%matbin=matbin';
matbin2=arrayfun(@(n)(sum(matbin(:,(n-1)*PRECI2+1:n*PRECI2),2))./(2^(PRECI-6)),1:size(ipop_new,2),'UniformOutput' , false);
TL_new_PSS=cell2mat(matbin2).*sgn;

  cost_new=zeros(size(   TL_new_PSS,1) ,1);
  
  for i=1:size(   TL_new_PSS,1) 
        cost_new(i,1) =  fun(TL_new_PSS(i,:));
  end
  
  costpss2=[cost_new,TL_new_PSS];
  costpss2=sortrows(costpss2,1);

  for m=1:numel(cost_new)
      [cs,ind]=max(Cost);
      if costpss2(m,1)<cs(1)
          %%% tavajoh kon ke inja leakage etelaat darim. (???)
          Cost(ind(1))= costpss2(m,1);
          Iop_PSS(ind(1),:)=costpss2(m,2:end);
      else
          break
      end
  end
        

  %  Cost = cost_new;
   %   Iop_PSS = TL_new_PSS;
    ipop_cost = [Iop_PSS,Cost];
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