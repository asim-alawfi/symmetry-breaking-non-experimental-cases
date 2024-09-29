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
%% Track symmetry breaking bifurcations for case TD=D=0.05
condproj=eye(xdim);
condS=eye(xdim);
pfxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.4,1)','stateproj',condS,'condprojmat',condS);
pofoldargs={'usercond',{pfxsym},'initcond',{pfxsym}};
[pffuncs_TDD,pfbranch_TDD,succ]=SetupPOEV1(funcs_audi,br_symmetry_wbifs(2),p_bif(2,2),...
    'contpar',[in.PR,in.df],'dir',in.PR,'step',0.1,'max_step',[in.PR,0.08; in.df,0.05; 0,0.05],...
     'print_residual_info',1,'use_tangent',0,'degree',6,'intervals',120,...
    pofoldargs{:});
%
%%
%pfbranch_TDD.parameter.max_step(end)=0.08;
figure(2)
clf 
hold on
pfbranch_TDD=br_contn(pffuncs_TDD,pfbranch_TDD,350);
%%
figure(2)
hold on
pfbranch_TDD=br_rvers(pfbranch_TDD);
pfbranch_TDD=br_contn(pffuncs_TDD,pfbranch_TDD,60);
%% Compute stability
[pfbranch_TDD_wbifs,nunst_tdd,bif_tdd,p_biftdd]=MonitorChange(pffuncs_TDD,pfbranch_TDD,...
    'range',2:length(pfbranch_TDD.point),'printlevel',1,'print_residual_info',0,...
    'min_iterations',5);
% parameters
rp_btdd=arrayfun(@(x)x.parameter(in.PR),pfbranch_TDD_wbifs.point);
df_btdd=arrayfun(@(x)x.parameter(in.df),pfbranch_TDD_wbifs.point);
%%
figure(23)
clf;
plot(rp_btdd(nunst_tdd==0),df_btdd(nunst_tdd==0),'k',...
     rp_btdd(nunst_tdd>=1),df_btdd(nunst_tdd>=1),'r.','LineWidth',2)
grid on
pfbranch_TDD=br_remove_extracolumns(pfbranch_TDD);
pfbranch_TDD_wbifs=br_remove_extracolumns(pfbranch_TDD_wbifs);
legend('symmetry-breaking','FontSize',16)
%%
save('symmetry_breaking_TD_and_D_try2.mat')
%% Plot symmetry-breaking + POs with touching theta %%%%%%%%%%%%%%
load('touching_theta_case_increased_TD_and_D_try2.mat')
%%
nunst_pct2=nunst_pc{2};
rp_theta2=arrayfun(@(x)x.parameter(in.PR),br_crossing_wbifs(2).point);
df_theta2=arrayfun(@(x)x.parameter(in.df),br_crossing_wbifs(2).point);

figure(23)
clf;
hold on
plot(rp_btdd(nunst_tdd==0),df_btdd(nunst_tdd==0),'k',...
     rp_btdd(nunst_tdd>=1),df_btdd(nunst_tdd>=1),'r.','LineWidth',2)
plot(rp_theta2(nunst_pct2==0),df_theta2(nunst_pct2==0),'m.',...
     rp_theta2(nunst_pct2>=1),df_theta2(nunst_pct2>=1),'r.','LineWidth',2)
plot(pr2_m(nunst2_mm==0),df2_m(nunst2_mm==0),'g.',...
    pr2_m(nunst2_mm>=1),df2_m(nunst2_mm>=1),'y.','LineWidth',2)
legend('symmetry-breaking','POs with touching \theta =0.5','','FontSize',16)
grid on
%
