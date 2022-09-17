function  mol_wt_mix  = mix_mol_wt(y)
%% Function to calculate the gas mixture molecular weight kg/mole
global prop_data;

mol_wt_mix = 0;
for i=1:length(y)-1
    mol_wt_mix = mol_wt_mix + y(i)* prop_data(i,1);
end

