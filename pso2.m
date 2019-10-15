function  [K3 , vK3, time1, solution, minimum ]=pso2(func,itr,subitr)
addpath 'my-optimizers\supp_functs'
evalc(['fun = @', func]);
[Nparam, bird_min_range, bird_max_range, solution, minimum] = feval(fun);
Nparam=numel(bird_min_range);
bird_min_range(bird_min_range==-Inf)=-10000;
bird_max_range(bird_min_range==Inf)=10000;

num_of_individuals=subitr;
max_iteration=itr;
format long


%% definitions 4 PSO
%MinMaxRange=[X_Min_PSS',X_Max_PSS'];
Bird_in_swarm=subitr;
Number_of_quality_in_Bird=Nparam;
availability_type='min';
velocity_clamping_factor=2;
cognitive_constant=2; %c1=individual learning rate (normally 2)
social_constant=2; %c2=social parameter (normally 2)
Min_Inertia_weight=0.4;
Max_Inertia_weight=0.9;

%% job
N=Bird_in_swarm*max_iteration;
q=0;
for i=1:Number_of_quality_in_Bird
    bird(:,i)=bird_min_range(i)+(bird_max_range(i)-bird_min_range(i))*rand(Bird_in_swarm,1);
end

Vmax=bird_max_range*velocity_clamping_factor;
Vmin=bird_min_range*velocity_clamping_factor;%-Vmax;

for i=1:Number_of_quality_in_Bird
    Velocity(:,i)=Vmin(i)+(Vmax(i)-Vmin(i))*rand(Bird_in_swarm,1);
end

for itr=1:max_iteration
   % fprintf('Completed  %d  %% ...', uint8(q*100/N ))
    tic
    for p=1:Bird_in_swarm
        parameter=bird(p,:,itr);
        availability(p,itr)=feval(fun,parameter);
        
        switch availability_type
            case 'min'
                format long;
                [pBest_availability,index]=min(availability(p,:));
                pBest=bird(p,:,index);
                
                if(p==1 && itr==1)
                    gBest=pBest;
                    gBest_availability=pBest_availability;
                elseif availability(p,itr)<gBest_availability
                    gBest_availability=availability(p,itr);
                    gBest=bird(p,:,itr);
                end
                
            case 'max'
                format long;
                [pBest_availability,index]=max(availability(p,:));
                pBest=bird(p,:,index);
                
                if(p==1 && itr==1)
                    gBest=pBest;
                    gBest_availability=pBest_availability;
                elseif availability(p,itr)>gBest_availability
                    gBest_availability=availability(p,itr);
                    gBest=bird(p,:,itr);
                end
                
            otherwise
                error('availability_type mismatch')
        end
        
        w(itr)=((max_iteration - itr)*(Max_Inertia_weight - Min_Inertia_weight))/(max_iteration-1) + Min_Inertia_weight;
        Velocity(p,:,(itr+1))=w(itr)*Velocity(p,:,itr) + social_constant*rand(1,Number_of_quality_in_Bird).*(gBest-bird(p,:,itr)) + cognitive_constant*rand(1,Number_of_quality_in_Bird).*(pBest-bird(p,:,itr));
        Velocity(p,:,(itr+1))=MinMaxCheck(Vmin, Vmax, Velocity(p,:,(itr+1)));
        
        bird(p,:,(itr+1))= bird(p,:,itr) + Velocity(p,:,(itr+1));
        bird(p,:,(itr+1))=MinMaxCheck(bird_min_range, bird_max_range, bird(p,:,(itr+1)));
        q=q+1;
    end
    
    clc;
    K3{itr}=gBest;
    vK3(itr)=fun(gBest);
    time1(itr)=toc;
end
optimised_parameters=gBest;