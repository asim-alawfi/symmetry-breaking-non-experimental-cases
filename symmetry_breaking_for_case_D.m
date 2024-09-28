%% Symmetry-breaking for the case when the delay is large $D=0.05$
clear;
    base=[pwd(),'\..\ddebiftool_snapshot_23October2022\'];
    base2=[pwd(),'\..\Supporting_function\'];
    addpath([base,'ddebiftool'],...
            [base,'ddebiftool_extra_psol'],...
            [base,'ddebiftool_utilities'],...
            [base,'ddebiftool_extra_rotsym'],...
            [base,'ddebiftool_extra_nmfm'],...
            [base,'ddebiftool_extra_symbolic'],...
            [base,'ddebiftool_coco']);
addpath([base2,'Supporting_functions'])
%%
 load('br_crossing_threshold_try2.mat')
%% Track symmetry breaking bifurcations for stst
%%
condproj=eye(xdim);
condS=eye(xdim);
pfxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.5,2)','stateproj',condS,'condprojmat',[1,0,0,0,0,0]);
pofoldargs={'usercond',{pfxsym},'initcond',{pfxsym}};
%%
bif_p=p_bif(1,3);
[pffuncs_D,pfbranch_D,succ]=SetupPOEV1(funcs_audi,br_symmetry_wbifs(3),p_bif(1,3),...
    'contpar',[in.PR,in.df],'dir',in.PR,'step',0.1,'max_step',[in.PR,0.1; in.df,0.05; 0,0.05],...
     'print_residual_info',1,'use_tangent',0,...
    pofoldargs{:});
%%
figure(2);
pfbranch_D=br_contn(pffuncs_D,pfbranch_D,179);
pfbranch_D.parameter.max_step(4)=0.01;
pfbranch_D=br_contn(pffuncs_D,pfbranch_D,30);
pfbranch_D.parameter.max_step(4)=0.1;
pfbranch_D=br_contn(pffuncs_D,pfbranch_D,218);
pfbranch_D=br_rvers(pfbranch_D);
pfbranch_D=br_contn(pffuncs_D,pfbranch_D,40);
%% Compute stability
[pfbranch_D_wbifs,nunst_sb,bif_b,p_bifb]=MonitorChange(pffuncs_D,pfbranch_D,...
    'range',2:length(pfbranch_D.point),'printlevel',1,'print_residual_info',0,...
    'min_iterations',5);
% parameters
rp_b=arrayfun(@(x)x.parameter(in.PR),pfbranch_D_wbifs.point);
df_b=arrayfun(@(x)x.parameter(in.df),pfbranch_D_wbifs.point);
%%
figure(33)
clf;
plot(rp_b(nunst_sb==0),df_b(nunst_sb==0),'r.',...
     rp_b(nunst_sb>=1),df_b(nunst_sb>=1),'k.','LineWidth',1.5)
grid on
%%
pfbranch_D=br_remove_extracolumns(pfbranch_D);
pfbranch_D_wbifs=br_remove_extracolumns(pfbranch_D_wbifs);
save('symmetry_breaking_D_try2.mat')
%% Starting tracking pitchfork(symmetry-breaking) from the other direction 
br_symmetry_wbifs(3).method.point.newton_max_iterations=5;
[pffuncs_D,pf_br,s]=SetupPOEV1(funcs_audi,br_symmetry_wbifs(3),p_bif(2,3),...
    'contpar',[in.PR,in.df],'dir',[],'step',0.01,'max_step',[in.PR,0.1; in.df,0.01; 0,0.005],...
     'print_residual_info',1,'use_tangent',0,...
    pofoldargs{:});
pf_br.method.point.newton_max_iterations=5;
new_point=pffuncs_D.get_comp(pf_br.point(1),'solution');
new_branch=br_symmetry_wbifs(3);
new_branch.point=new_point;
[pffuncs_D,pf_br,suc]=SetupPOEV1(funcs_audi,new_branch,1,...
    'contpar',[in.PR,in.df],'dir',in.PR,'step',0.1,'max_step',[in.PR,0.1; in.df,0.05; 0,0.05],...
     'print_residual_info',1,'use_tangent',0,...
    pofoldargs{:});
figure(2)
pf_br=br_contn(pffuncs_D,pf_br,26);
pf_br=br_rvers(pf_br);
pf_br=br_contn(pffuncs_D,pf_br,200);
%%
[pf_br_wbifs,uns,bif_b2,p_bifb2]=MonitorChange(pffuncs_D,pf_br,...
    'range',2:length(pf_br.point),'printlevel',1,'print_residual_info',0,...
    'min_iterations',5);
rp_b2=arrayfun(@(x)x.parameter(in.PR),pf_br_wbifs.point);
df_b2=arrayfun(@(x)x.parameter(in.df),pf_br_wbifs.point);
%% Plot the symmetry-breaking + the touching branches
%load('br_crossing_threshold_try2.mat')
%%
nunst_pc3=nunst_pc{3};
rp_bt2=arrayfun(@(x)x.parameter(in.PR),br_crossing_wbifs(3).point);
df_bt2=arrayfun(@(x)x.parameter(in.df),br_crossing_wbifs(3).point);
figure(34)
clf
hold on 
plot(rp_b(nunst_sb==0),df_b(nunst_sb==0),'r.',...
     rp_b(nunst_sb>=1),df_b(nunst_sb>=1),'k.','LineWidth',1.5)
plot(rp_b2(uns==0),df_b2(uns==0),'r.',...
     rp_b2(uns>=1),df_b2(uns>=1),'k.','LineWidth',1.5)
plot(rp_bt2(nunst_pc3==0),df_bt2(nunst_pc3==0),'m--','LineWidth',2)
plot(rp_bt2(nunst_pc3>=1),df_bt2(nunst_pc3>=1),'.','Color',[0.7 0.7 0.7],'LineWidth',2)
grid on
%%
pf_br=br_remove_extracolumns(pf_br);
pf_br_wbifs=br_remove_extracolumns(pf_br_wbifs);
%%
save('symmetry_breaking_D_b2_try2.mat')
%%
load('touching_theta_case_increased_D_try2.mat')
nunst_pc3=nunst_pc{3};
rp_bt2=arrayfun(@(x)x.parameter(in.PR),br_crossing_wbifs(3).point);
df_bt2=arrayfun(@(x)x.parameter(in.df),br_crossing_wbifs(3).point);
%%
figure(34)
clf
hold on 
plot(rp_b(nunst_sb==0),df_b(nunst_sb==0),'r.',...
     rp_b(nunst_sb>=1),df_b(nunst_sb>=1),'k.','LineWidth',1.5)
plot(rp_b2(uns==0),df_b2(uns==0),'r.',...
     rp_b2(uns>=1),df_b2(uns>=1),'k.','LineWidth',1.5)
plot(rp_bt2(nunst_pc3==0),df_bt2(nunst_pc3==0),'m--','LineWidth',2)
plot(rp_bt2(nunst_pc3>=1),df_bt2(nunst_pc3>=1),'.','Color',[0.7 0.7 0.7],'LineWidth',2)
plot(pr3_m(nunst3_mm==0),df3_m(nunst3_mm==0),'k.',...
    pr3_m(nunst3_mm>=1),df3_m(nunst3_mm>=1),'y.','LineWidth',2)
grid on