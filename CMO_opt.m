%% starting CMO code:
function [K3 , vK3 , time1, solution, minimum ]=CMO_opt(func,itr,subitr)
%% general params
addpath 'my-optimizers\supp_functs'
evalc(['fun = @', func]);
global termination_of_criteria %%% (TODO)
global conditions_value
global conditions_ordered_by_prefered_priorities
global solutions
global costs
global conditions_value_gama
global maximum_conflict_log
global defaultprognum
global max_num_solutions
global alpha
tic

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

%% predefine spec
THRESH0=1.2; %%% (TODO)
num_of_individuals=subitr;
max_data=ceil(num_of_individuals*THRESH0);
all_data=zeros(num_of_individuals,Nparam); %? 
all_cost=zeros(num_of_individuals,1); %?
scores=ones(num_of_individuals,1);
all_scores=ones(max_data,num_of_individuals);
format long

%% now, define the main specs
%%% (TODO) These have to be set:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conditions_value=[]; %(TODO) fill in it asap (each condition has specific value which gives that priority)
conflict_flag=0; last_conflicted_condition=0; conflicted_condition=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

programs_options={}; %(TODO) set its default
conditions_options={}; %(TODO) set its default
conditions_options{1}.func=func;

each_programs_default_variables={}; %%% (TODO) set them , set min of them set max of them
% each condition should have series of parameters and its favorite prgram
program_num_elected_4_run=defaultprognum;

% %%% neeeded params to put into options 1 unwanted (todo)
% programs_options{1}.xmin=xmin;programs_options{1}.xmmax=xmax;

% preferedprognum ??? (TODO)
% any of them should not have their default options removed but saved, even
% after passing the update_prog_options

%%%%%%
programs_name=dir('my-optimizers\CMO_depends\programs_archive\*.m');
programs_name = cellfun(@(x) x(1:end), {programs_name.name}, 'uniformoutput', false);
programs_name=regexp(programs_name,'.m','split');
programs_name=arrayfun(@(m)programs_name{m}{1},1:numel(programs_name),'UniformOutput', false);
programs_len=numel(programs_name);
addpath('my-optimizers\CMO_depends\programs_archive\');
%%%%%%
conditions_name=dir('my-optimizers\CMO_depends\if_then_archive\*.m');
conditions_name = cellfun(@(x) x(1:end), {conditions_name.name}, 'uniformoutput', false);
conditions_name=regexp(conditions_name,'.m','split');
conditions_name=arrayfun(@(m)conditions_name{m}{1},1:numel(conditions_name),'UniformOutput', false);
conditions_len=numel(conditions_name);
addpath('my-optimizers\CMO_depends\if_then_archive\');
%%%%%%
reactions_name = dir('my-optimizers\CMO_depends\if_then_archive\*.m');
reactions_name = cellfun(@(x) x(1:end), {reactions_name.name}, 'uniformoutput', false);
reactions_name = regexp(reactions_name,'.m','split');
reactions_name = arrayfun(@(m)reactions_name{m}{1},1:numel(reactions_name),'UniformOutput', false);
addpath('my-optimizers\CMO_depends\if_then_archive\');

%% now call options functions to fill in prerequisites
addpath('my-optimizers\CMO_depends\');
[programs_options,conditions_options]=call4options(programs_options,conditions_options,func); %(TODO)

%% randomize initial solutions
for i=1:num_of_individuals
    rnd=rand(1,numel(xmin(1,:)));
    ipop_new(i,1:Nparam)=  rnd.*(xmin(1,:)-xmin(1,:))+(1-rnd).*xmin(1,:);
    cost_new(i,1)=fun( ipop_new(i,1:Nparam));
    update_alldata( ipop_new(i,1:Nparam),cost_new(i,1) ,i)
end
ipop_cost = [ipop_new,cost_new];
ipop_cost_sort = sortrows(ipop_cost,Nparam+1);
solutions = ipop_cost_sort(:,1:Nparam);
costs = ipop_cost_sort(:,Nparam+1);


%% ask for all existing programs  from program_archive
%programs= load_all_programs( from_archive);
% ask for all existing rules from rule_archive
Iter=0;
while and(~termination_of_criteria, itr>Iter)
    Iter=Iter+1;
    m=0;
    conflict_happened4one_iter=zeros(1,conditions_len);
    conflict_intensity4one_iter=zeros(1,conditions_len);
    %     conflict_intensity(1:end-maximum_conflict_log,:)=[]; % (todo) MAYBE DEPRECATED % for restrincting its size
    
    while ~conflict_flag
        m=m+1;
        %% conditioning phase
        [solutions, costs, programs_options{program_num_elected_4_run}] = feval(programs_name{program_num_elected_4_run}, solutions, costs, func, programs_options{program_num_elected_4_run});
        j=0;
        
        % check whether or not solution count exceeded than number it deserves
        if size(solutions,1)>max_num_solutions
            %             solutions(1:end-max_num_solutions,:)=[];
            %             costs(1:end-max_num_solutions,:)=[];
            %              solutions(1+max_num_solutions:end,:)=[];
            %             costs(1+max_num_solutions:end,:)=[];
            [~,inn]=sort(costs);solutions=solutions(inn(1:max_num_solutions),:);costs=costs(inn(1:max_num_solutions));
            
            
        end
        
        %% autonomic nervous system , unCMO conflict verification
        % First of all, associate each condition to its related programs
        
        %%
        if conflicted_condition~=0, last_conflicted_condition=conflicted_condition; end
        conflicted_condition=0;
        for i=conditions_ordered_by_prefered_priorities % as mohem be kam ahamiat
            j=j+1;
            if sum(conditions_options{i}.validprograms==program_num_elected_4_run) > 0
                [conflict_happened4one_iter(j), conflict_intensity4one_iter(i), programs_options{program_num_elected_4_run}, conditions_options{i}]=feval(conditions_name{i}, programs_options{program_num_elected_4_run}, conditions_options{i});
                
                % conflict occured
                if conflict_happened4one_iter(j)==1,
                    conflicted_condition=i;
                    conditions_options{i}.intensity_before_reaction=[conditions_options{i}.intensity_before_reaction,conflict_intensity4one_iter(i)];
                    break;
                end
            end
        end
                
        
        
        %  if conflicted_condition==0, termination_of_criteria=1; % (todo) Set other terminations of criteria
        %         conflict_happened=[conflict_happened;conflict_happened4one_iter(conditions_ordered_by_prefered_priorities)]; % (todo) MAYBE DEPRECATED
        %         conflict_intensity=[conflict_intensity;conflict_intensity4one_iter]; % (todo) MAYBE DEPRECATED
        if conflicted_condition >0 ,conflict_flag=1; end
    end
    %% now update/modify the selected program's options
    programs_options= feval(reactions_name{conflicted_condition},programs_options, conditions_options{conflicted_condition});
    
    % each condition should have series of parameters and its favorite prgram
    program_num_elected_4_run1=conditions_options{conflicted_condition}.preferedprognum;
    if program_num_elected_4_run1>0 % for those cases of conditions without affecting program parameters
        program_num_elected_4_run=program_num_elected_4_run1;
    end
    % then tune the very conditions' floating parameters (???) (TODO)
    
    G_best_value=costs(end);G_best=solutions(end,:);
    Nparam=size(solutions,2);
    ipop_cost = [solutions,costs];
    ipop_cost_sort=sortrows(ipop_cost,Nparam+1);
    if ipop_cost_sort(1,Nparam+1) < G_best_value
        G_best = ipop_cost_sort(1,1:Nparam);
        G_best_value = ipop_cost_sort(1,Nparam+1);
    end
    %  Tr = [G_best,G_best_value];
    K3{Iter} = G_best;
    vK3(Iter) = G_best_value;
    time1(Iter)=toc;
    conflict_flag=0;
    
end

%time1=toc;
time1

