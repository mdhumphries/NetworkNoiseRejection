% gather all results into a Table...
clear all; close all;

fnames = dir('../Results/');

nF = numel(fnames);

netCtr = 0;
for iF = 1:nF
    
    if any(strfind(fnames(iF).name,'Rejected'))
        netCtr = netCtr + 1;
        result(netCtr).NetworkName = fnames(iF).name(10:end-4); % strip out 'Rejected' and .mat
        load(['../Results/' fnames(iF).name]);
        result(netCtr).Network_Size = numel(Data.ixRetain); 
        % network density
        % keyboard
        result(netCtr).Network_Density = sum(sum(Data.A>0)) / (result(netCtr).Network_Size.^2 - result(netCtr).Network_Size);  % undirected, no self-loops 
        % network type
        if all(Data.A(Data.A>0) == 1)
            result(netCtr).LinkType = 'binary';
        elseif ~any(rem(Data.A(:),1))  % then is integers for all weights
            result(netCtr).LinkType = 'integer';
        else
            result(netCtr).LinkType = 'real';
        end
        % results of rejection
        result(netCtr).Signal_Size_Sparse = numel(Data.ixSignal_Final);   
        result(netCtr).Signal_Size_Full = numel(Control.ixSignal_Final); 
        result(netCtr).SparseWCM_Dn = Data.Dn;
        result(netCtr).FullWCM_Dn = Control.Dn;
        result(netCtr).SparseWCM_PosDn = Data.PosDn;
        result(netCtr).FullWCM_PosDn = Control.PosDn;
        result(netCtr).Signal_Components_Sparse = numel(Data.SignalComp_sizes);
        % result(netCtr).Signal_Components_Full = numel(Control.SignalComp_sizes);
        
        % negative eigenvalues
        result(netCtr).SparseWCM_NegDn = Data.Dneg;
        result(netCtr).FullWCM_NegDn = Control.Dneg;
        
        % basic properties
        property(netCtr).Name = fnames(iF).name(10:end-4);
        property(netCtr).N = result(netCtr).Network_Size;
        property(netCtr).m = sum(sum(Data.A>0));
        property(netCtr).density = result(netCtr).Network_Density;
        property(netCtr).LinkType = result(netCtr).LinkType;
    end
end

Network_Rejection_Table = struct2table(result);
save('../Results/Network_Rejection_Table','Network_Rejection_Table');

Network_Property_Table = struct2table(property);
save('../Results/Network_Property_Table','Network_Property_Table');

