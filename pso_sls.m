function  [K3 , vK3, time1, solution, minimum ]=pso_sls(func,itr,subitr)
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

%% definitions 4 PSO
%MinMaxRange=[X_Min_PSS',X_Max_PSS'];
Bird_in_swarm=num_of_individuals;
Number_of_quality_in_Bird=Nparam;
%availability_type='min';
velocity_clamping_factor=2;
cognitive_constant=2; %c1=individual learning rate (normally 2)
social_constant=2; %c2=social parameter (normally 2)
Min_Inertia_weight=0.4;
Max_Inertia_weight=0.9;

%% initial
N=Bird_in_swarm*max_iteration;
q=0;
for i=1:Bird_in_swarm
    bird(i,:)=xmin+(xmax-xmin).*rand(1,Number_of_quality_in_Bird);
    parameter=bird(i,:);
    availability(i)=feval(fun,parameter);
    update_alldata(parameter,availability(i) ,i)
end
pBest = bird;pBestAvail=availability;

[~,indx]=min(availability);
gBest = bird(indx,:);gBestAvail=availability(indx);

Vmax=xmax*velocity_clamping_factor;
Vmin=xmin*velocity_clamping_factor;%-Vmax;

for i=1:Number_of_quality_in_Bird
    Velocity(:,i)=Vmin(i)+(Vmax(i)-Vmin(i))*rand(Bird_in_swarm,1);
end




for itr=1:max_iteration
   % fprintf('Completed  %d  %% ...', uint8(q*100/N ))
    tic
    
    %% job 4 pso
    for p=1:Bird_in_swarm
        
        
        % update solution
        w(itr)=((max_iteration - itr)*(Max_Inertia_weight - Min_Inertia_weight))/(max_iteration-1) + Min_Inertia_weight;
        Velocity(p,:)=w(itr)*Velocity(p,:) + social_constant*rand(1,Number_of_quality_in_Bird).*(gBest-bird(p,:)) + cognitive_constant*rand(1,Number_of_quality_in_Bird).*(pBest(p,:)-bird(p,:));
        Velocity(p,:)=MinMaxCheck(Vmin, Vmax, Velocity(p,:));
        
        lbird=bird(p,:); bird(p,:)= lbird + Velocity(p,:);
        bird(p,:)=MinMaxCheck(xmin, xmax, bird(p,:));
      %  q=q+1;
        
        
        % update cost
        lavail =availability(p); availability(p)= fun(bird(p,:));
        
        % update pBest, gBest
        if availability(p)<pBestAvail(p)
            pBest(p,:) = bird(p,:);
            pBestAvail(p)=availability(p);
            update_alldata(bird(p,:), availability(p) ,p)
        else
            % baresh gardun!
            update_alldata(bird(p,:), availability(p) ,0)
            bird(p,:)=lbird;
            availability(p)=lavail; 
        end    
        if availability(p)<gBestAvail
            gBest= bird(p,:);
            gBestAvail=availability(p);
        end
        
 
    end
    
    
    
    
%     clc; %???
%     K3{itr}=gBest; %???
%     vK3(itr)=fun(gBest); %???
%     time1(itr)=toc; %???

% %% get data ready for sls part
% Cost=availability;
% Iop_PSS=bird;

    %% job 4 sls
    for i = 1:num_of_individuals

        LS_new2=update_linear();
 

        costnew2 = fun(LS_new2);

        if costnew2 < availability(sel1)
            scores(sel1)=scores(sel1)/THRESH4SCORE;
            all_scores(:,sel1)=1;
          if sel2~=0
              all_scores(sel2,sel1)=all_scores(sel2,sel1)/THRESH4SCORE; 
          end
            
            bird(sel1,1:Nparam) =LS_new2; %MBC
            availability(sel1) = costnew2; %MBC
            update_alldata(LS_new2,costnew2 ,sel1)
           %%%%%%%%%55
            p=sel1;
            if availability(p)<pBestAvail(p)
                pBest(p,:) = bird(p,:);
                pBestAvail(p)=availability(p);
            end    
            if availability(p)<gBestAvail
                gBest= bird(p,:);
                 gBestAvail=availability(p);
            end
           
           %%%%%%%%%%%

       else
            update_alldata( LS_new2, costnew2 ,0)
            scores(sel1)=scores(sel1)*THRESH4SCORE;
             if sel2~=0
            all_scores(sel2,sel1)=all_scores(sel2,sel1)*THRESH4SCORE;
             end
       end
    end

%     %% send result of SLS to PSO %MBC
%     bird(sel1,:,itr)=LS_new2;
%     availability(sel1,itr)=costnew2;
            
%     ipop_cost = [bird,cost_new];
%     ipop_cost_sort=sortrows(ipop_cost,Nparam+1);
%     if ipop_cost_sort(1,Nparam+1) < G_best_value
%         G_best = ipop_cost_sort(1,1:Nparam);
%         G_best_value = ipop_cost_sort(1,Nparam+1);
%     end
    %Tr = [G_best,G_best_value];
        K3{itr} = gBest;
vK3(itr) = gBestAvail;
 time1 (itr)=toc;





end


 optimised_parameters=gBest;