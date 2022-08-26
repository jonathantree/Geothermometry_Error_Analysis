%% Jonathan Tree

%This script conducts simple olivine removal by finding the median
%forsterite measured within a sample and subtracting small steps from the
%whole rock composition to obtain the parental magma in equilibrium with
%the highest fosterite measured in that sample
clear all; clc;

% sample_name='NIH-F-5A';
% 
workbookFile = 'D18_1_ol.csv';
% sheetName='sheet1';
% startRow=2;
% endRow=27;
% 
% [SiO2_Ol, FeO_Ol, MgO_Ol, CaO_Ol, NiO_Ol, MnO_Ol, Total, Fo_Ol] = ...
%     import_ol_csv(workbookFile);
%%
% [SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O3 K2O P2O5 Total]
WR = [45.09	3.59	16.62	(0.1*14.33) (0.9*14.33)	0.13	4.49	10.84	3.11	0.68	1.10	99.97];
%%
median_Fo = median(Fo_Ol);
mean_Fo   = mean(Fo_Ol);
mode_Fo   = mode(Fo_Ol);
n         = length(Fo_Ol);


for i= 1:n
    tmp(i) = abs(median_Fo-Fo_Ol(i));
    tmp=tmp';
    [x, idx] = min(tmp); %index of closest value
    closest_Fo = Fo_Ol(idx); %closest Fo
end

display (['The Olivine Forsterite being taken out is ' num2str(closest_Fo)])

Model_OL_SiO2 = SiO2_Ol(idx);
Model_OL_FeO  = FeO_Ol(idx);
Model_OL_NiO  = NiO_Ol(idx);
Model_OL_MnO  = MnO_Ol(idx);
Model_OL_MgO  = MgO_Ol(idx);
Model_OL_CaO  = CaO_Ol(idx);

N_WR        = length(WR);
Model_Magma(1, 1:N_WR) = WR;   

% [SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O3 K2O P2O5 Total]

OL_Model_Array = [Model_OL_SiO2 0 0 0 Model_OL_FeO Model_OL_MnO Model_OL_MgO ...
                    Model_OL_CaO 0 0 0 0];

%%

MAX_Fo    = max(Fo_Ol);
display(['The Max Fo to equilibrate with is ' num2str(MAX_Fo)])


Kd1 = 0.315; 


Mg_Fe_Ol   = (1:0.1:1000)';
Mg_NUM_ol  = 100*(Mg_Fe_Ol./(Mg_Fe_Ol+1))';

Mg_Fe_liq_1  = (Mg_Fe_Ol*Kd1)';
Mg_NUM_liq_1 = 100*(Mg_Fe_liq_1./(Mg_Fe_liq_1+1))';
    

for o= 1:length(Mg_NUM_ol)
    Fo_tmp(o) = abs(MAX_Fo-Mg_NUM_ol(o));
    [Fo, idx_Max_Fo] = min(Fo_tmp); %index of closest value
    Target_MG_Num = Mg_NUM_liq_1(idx_Max_Fo); %closest Fo
end

%%
LIMIT = 0.01;
flag = 0;
d_Mg(2) = 5;
j = 2;

while (flag==0) 
    Model_Magma(j, 1:N_WR) = Model_Magma(j-1, 1:N_WR) - 0.0001*OL_Model_Array;
    TOTAL(j) = sum(Model_Magma(j, 1:N_WR-1));
    Model_Magma(j, N_WR) = TOTAL(j);
    Model_Magma(j, 1:N_WR) = 100*(Model_Magma(j, 1:N_WR)./TOTAL(j));
    Model_Magma(j, N_WR) = sum(Model_Magma(j, 1:N_WR-1));
    
        Magma_Mg(j) = Model_Magma(j, 7)/40.31;
        Magma_Fe(j) = Model_Magma(j, 5)/71.85;

        Magma_Mg_Num(j) = (Magma_Mg(j)./(Magma_Mg(j)+Magma_Fe(j)))*100;
        
            
        d_Mg(j) = abs(Magma_Mg_Num(j)-Target_MG_Num);
        
        if d_Mg(j) < LIMIT
            flag = 1;
        elseif d_Mg(j) > LIMIT
            j=j+1;
        end 
      
end  
[Kd_315] = Model_Magma(j, 1:N_WR);
Kd_315 = Kd_315';

%% Find the parental magmas in equilibrium with Kd = 0.345 and 0.375
Kd2 = 0.345;
Mg_Fe_liq_2  = (Mg_Fe_Ol*Kd2)';
Mg_NUM_liq_2 = 100*(Mg_Fe_liq_2./(Mg_Fe_liq_2+1))';

for oo= 1:length(Mg_NUM_ol)
    Fo_tmp(oo) = abs(MAX_Fo-Mg_NUM_ol(oo));
    [Fo, idx_Max_Fo] = min(Fo_tmp); %index of closest value
    Target_MG_Num2 = Mg_NUM_liq_2(idx_Max_Fo); %closest Fo
end

for k= 1:length(Magma_Mg_Num)
    Mg_Num2_tmp(k) = abs(Target_MG_Num2-Magma_Mg_Num(k));
    [Mg_num2, idx_Mg_num2] = min(Mg_Num2_tmp); %index of closest value
    [Kd_345] = Model_Magma(idx_Mg_num2, 1:N_WR); %closest Fo
    Kd_345 = Kd_345';
end

Kd3 = 0.375;
Mg_Fe_liq_3  = (Mg_Fe_Ol*Kd3)';
Mg_NUM_liq_3 = 100*(Mg_Fe_liq_3./(Mg_Fe_liq_3+1))';

for ooo= 1:length(Mg_NUM_ol)
    Fo_tmp(ooo) = abs(MAX_Fo-Mg_NUM_ol(ooo));
    [Fo, idx_Max_Fo] = min(Fo_tmp); %index of closest value
    Target_MG_Num3 = Mg_NUM_liq_3(idx_Max_Fo); %closest Fo
end

for kk= 1:length(Magma_Mg_Num)
    Mg_Num3_tmp(kk) = abs(Target_MG_Num3-Magma_Mg_Num(kk));
    [Mg_num3, idx_Mg_num3] = min(Mg_Num3_tmp); %index of closest value
    [Kd_375] = Model_Magma(idx_Mg_num3, 1:N_WR); %closest Fo
    Kd_375 = Kd_375';
end


%% Display the Results

display (['The target Mg # for equilibrium with Kd=0.315 is ' num2str(Target_MG_Num)]);
display (['The final magma Mg # for Kd=0.315 is ' num2str(Magma_Mg_Num(j))]);

display (['The target Mg # for equilibrium with Kd=0.345 is ' num2str(Target_MG_Num2)]);
display (['The final magma Mg # for Kd=0.345 is ' num2str(Magma_Mg_Num(idx_Mg_num2))]);


display (['The target Mg # for equilibrium with Kd=0.375 is ' num2str(Target_MG_Num3)]);
display (['The final magma Mg # for Kd=0.375 is ' num2str(Magma_Mg_Num(idx_Mg_num3))]);


Oxides={'SiO2'; 'TiO2'; 'Al2O3'; 'Fe2O3'; 'FeO'; 'MnO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'P2O5'; 'Total'};
Parental_Magamas=table(Kd_315, Kd_345, Kd_375,  'RowNames', Oxides)

writetable (Parental_Magamas, 'Ol_Removal_PM.xlsx', 'WriteRowNames', true, 'Sheet', sample_name);








































