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
%% Track symmetry breaking bifurcations for 
condproj=eye(xdim);
condS=eye(xdim);
pfxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.4,1)','stateproj',condS,'condprojmat',condS);
pofoldargs={'usercond',{pfxsym},'initcond',{pfxsym}};
[pffuncs_TD,pfbranch_TD,succ]=SetupPOEV1(funcs_audi,br_symmetry_wbifs(1),p_bif(2,1),...
    'contpar',[in.PR,in.df],'dir',in.PR,'step',0.1,'max_step',[in.PR,0.05; in.df,0.05; 0,0.05],...
     'print_residual_info',1,'use_tangent',0,'degree',6,'intervals',120,...
    pofoldargs{:});
%
pfbranch_TD.parameter.max_bound(3)=28;
%%
figure(2)
hold on
%clf;
pfbranch_TD=br_contn(pffuncs_TD,pfbranch_TD,450);
%%
pfbranch_TD.parameter.max_step(4)=0.01;
pfbranch_TD.parameter.max_step(5)=0.01;
figure(2)
hold on
pfbranch_TD=br_contn(pffuncs_TD,pfbranch_TD,30);
%%
pfbranch_TD.parameter.max_step(4)=0.05;
pfbranch_TD.parameter.max_step(5)=0.05;
figure(2)
hold on
pfbranch_TD=br_contn(pffuncs_TD,pfbranch_TD,500);
%%
pfbranch_TD=br_rvers(pfbranch_TD);
pfbranch_TD=br_contn(pffuncs_TD,pfbranch_TD,60);
%% Compute stability
[pfbranch_TD_wbifs,nunst_td,bif_td,p_biftd]=MonitorChange(pffuncs_TD,pfbranch_TD,...
    'range',2:length(pfbranch_TD.point),'printlevel',1,'print_residual_info',0,...
    'min_iterations',5);
% parameters
rp_btd=arrayfun(@(x)x.parameter(in.PR),pfbranch_TD_wbifs.point);
df_btd=arrayfun(@(x)x.parameter(in.df),pfbranch_TD_wbifs.point);
%%
figure(13)
clf;
plot(rp_btd(nunst_td==0),df_btd(nunst_td==0),'r.',...
     rp_btd(nunst_td>=1),df_btd(nunst_td>=1),'k.','LineWidth',1.5)
grid on
%
pfbranch_TD=br_remove_extracolumns(pfbranch_TD);
pfbranch_TD_wbifs=br_remove_extracolumns(pfbranch_TD_wbifs);
save('symmetry_breaking_TD_try2.mat')
%% %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
