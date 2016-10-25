%% params
%static
cc=[100, 20];


%looping
fraction_periphery_grid=[0.2 0.5 0.8 1.0];
gg_cc_grid=[0.1 0.3 0.5 0.7 0.9];
cp_gg_grid=[0.1 0.3 0.5 0.7 0.9];

%% variables
core_perf=table;


%% set loops for params
for count_peri=1:numel(fraction_periphery_grid)
    fraction_periphery=fraction_periphery_grid(count_peri);
    for count_gg_cc=1:numel(gg_cc_grid)
        gg_cc=gg_cc_grid(count_gg_cc);
        mixing_index=gg_cc;
        for count_cp_gg=1:numel(cp_gg_grid)
            cp_gg=cp_gg_grid(count_cp_gg);
            noise_index=cp_gg;
            mixing=mixing_index * cc;
            noise=noise_index * mixing; 
            pp=noise/2;
            %generate planted noise network
            options.weight_dist=struct('ingroup', [100,20], 'outgroup', [1,0.02],...
                'periphery_periphery', [1,0.02], 'periphery_core', [10,2]);
            network=test_noise_rejection_planted_noise(50,2, 'low',1.0, options);
            % unit analyses
            
            signal_nodes=reject_the_noise(network.adjacency,'temp_planted_noise');
            true_members=network.membership;
            true_members(true_members >=0)=1;
            true_members(true_members<0)=0;
            performance=sum(true_members==signal_nodes)/numel(true_members);
            
            temp_table=table(fraction_periphery, mixing_index, noise_index, performance);
            core_perf=[core_perf; temp_table];
            
        end
    end

end
%% visualise loop analyses

