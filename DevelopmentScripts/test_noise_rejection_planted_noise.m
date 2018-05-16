function network=test_noise_rejection_planted_noise(nodes_group,num_groups, mix,frac_periphery,options)
%%% generate weighted stochastic block model network for testing noise
%%% rejection
%%%
%%% Inputs:
%%%     nodes_group: no. of nodes /group
%%%     num_groups: no. of communities
%%%     mix: 'low', 'medium', 'high'
%%%     periphery: fraction of nodes relative to core that are peripheral
%%%     options: structure of options
%%%           .blnViz: (0,1): plot resulting weight matrix (1 is default)
%%%           .weight_dist: struct with following fields (this overwrites
%%%           settings in "mix")
%%%                  .ingroup
%%%                  .outgroup
%%%                  .periphery_periphery
%%%                  .periphery_core
%%%            
%%% Outputs:
%%%     network: struct with following fields
%%%                  adjacency: adjacency matrix
%%%                  membership: community index [0:num_groups-1], peripheral nodes represented by -1 
%%% EXAMPLE: 
%%%      test_noise_rejection_planted_noise(50,2, 'low',0.2, options) 
%%%      where 
%%%      options.weight_dist=struct('ingroup', [100,20], 'outgroup', [1,0.02],...
%%%        'periphery_periphery', [1,0.02], 'periphery_core', [10,2]);
%%%      produces
%%%      network with 2 communities, each with 50 nodes and low mixing
%%%      between the communities and 0.2 * (50 * 2) peripheral noise nodes
%%%     
%%%      test_noise_rejection_planted_noise([60 40 82],3,'low', 0.3, options)
%%%      produces network with  3 communities of sizes 60, 40 and 82
% Abhinav Singh & Mark Humphries

addpath('./WSBM_v1.2/'); 
addpath('./WSBM_v1.2/analysis tools/')
% if nargin==0,
%     
%     nodes_group=50;
%     num_groups=2;
%     mix='low';
%     frac_periphery=0.5;
%     
%     weight_dist=struct('ingroup', [100,20], 'outgroup', [1,0.02],...
%         'periphery_periphery', [1,0.02], 'periphery_core', [10,2]);
%      
% end
%%
% R = [1,2,3,3;
%     2,4,3,3;
%     3,3,1,2;
%     3,3,2,4];
% theta_w = [5,1;20,1; 0,1; 100,1];
% theta_e = [1; 1; 1; 1];
% group_sizes = [25;25;25;25];
% [~,True_Model] = generateEdges('Normal','Bernoulli',...
%     R,theta_w,theta_e,group_sizes);
%% neurons in periphery - core -periphery- core

%num_groups=2;

%nodes_group=floor(N/num_groups);
if numel(nodes_group)>1
    num_periphery=floor(frac_periphery * sum(nodes_group));
    group_sizes=[nodes_group'; num_periphery];
    num_groups=numel(nodes_group);
else
    num_periphery=floor(frac_periphery* (num_groups * nodes_group)); 

    %num_core=nodes_group-num_periphery;
    group_sizes=[nodes_group*ones(num_groups,1); num_periphery];
end

num_nodes=sum(group_sizes);

%% edge existence probability
theta_e=ones(4,1);

%% weight distribution probability
%theta_w=zeros(numel(unique(R)),2);
if isfield(options, 'weight_dist')
    weight_dist=options.weight_dist;
    core_core=weight_dist.ingroup;
    group_group=weight_dist.outgroup;
    periphery_periphery=weight_dist.periphery_periphery;
    periphery_core=weight_dist.periphery_core;
    
else


    if strcmp(mix,'low')
        core_core=[100, 20];
        periphery_core=core_core/10;
        periphery_periphery=core_core/100;
        group_group=core_core/100;
        
    elseif strcmp(mix, 'medium')
        core_core=[100, 20];
        periphery_core=core_core/10;
        periphery_periphery=core_core/10;
        group_group=periphery_periphery;
        
    elseif strcmp(mix, 'high')
        core_core=[100, 20];
        periphery_core=core_core/1.5;
        periphery_periphery=core_core/2;a
        group_group=core_core/2;
        
    end
end

if ~isfield(options, 'blnViz')
    options.blnViz = 1;
end

% theta_w=[periphery_periphery;periphery_core;
%         group_group;core_core];

theta_w=[core_core; group_group; periphery_periphery; periphery_core];
%% weight family label for each connection between blocks
% R = [1,2,3,3;
%     2,4,3,3;
%     3,3,1,2;
%     3,3,2,4];

temp_R=2*ones(num_groups) -eye(num_groups);
R=zeros(num_groups+1);
R(1:end-1,1:end-1)=temp_R;
R(end,end)=3;
R(end,1:end-1)=4;
R(1:end-1,end)=4;
%% generate network by edge list
[~,True_Model] = generateEdges('Normal','Bernoulli',...
    R,theta_w,theta_e,group_sizes);
%% plot edge matrix 
if options.blnViz
plotWSBM(True_Model);
end

%% convert edge matrix into symmetric adjacency matrix
adj=zeros(num_nodes);
edges=True_Model.Data.Raw_Data;
ind=sub2ind(size(adj), edges(:,1), edges(:,2));
adj(ind)=edges(:,3);
adj_tiu=triu(adj);
adj=adj_tiu + transpose(adj_tiu) - diag(diag(adj_tiu));

network.adjacency=adj;

%% group membership
membership=zeros(num_nodes,1);
for i=0:num_groups-1
    membership(i*nodes_group+1:(i+1)*nodes_group,1)=i;
end
membership(num_groups*nodes_group+1:end)=-1;
network.membership=membership;










