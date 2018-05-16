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
        result(netCtr).Signal_Size_WCM = numel(Data.ixSignal_Final);        
        result(netCtr).WCM_Dn = Data.PosDn;
        result(netCtr).WCM_RejectionDn = Data.Dn;
        result(netCtr).Config_Dn = Control.PosDn;
        result(netCtr).Config_RejectionDn = Control.Dn;
        result(netCtr).Signal_Components = numel(Data.SignalComp_sizes);
    end
end

Network_Rejection_Table = struct2table(result);
save('../Results/Network_Rejection_Table','Network_Rejection_Table');