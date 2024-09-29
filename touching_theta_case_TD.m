%% We compute a branches of periodic solutions (symmetric and non-symmetric) by varying df, and we pick
% a non-symmetric solution with touching a threshold, then we continue
% the computation tracking a branch of non-symmetric solutions touching the
% threshold.
%%
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
%% One-parameter continuation with fixing PR=12, and varying df.
parbds={'min_bound',[in.PR,1;in.df,0.005],'max_bound',[in.PR,34.7; in.df,1],...
    'max_step',[in.PR,0.01; in.df,0.01; 0,0.005],'print_residual_info',1};
[~,it]=min(abs(rp1_x-7));
[br1_df,sucp]=ChangeBranchParameters(funcs_ss,br_symmetry_wbifs(1),it,...
      'degree',6,'intervals',120,'contpar',in.df,'step',1e-3,'correc',true,parbds{:},...
      'extra_condition',true,'phase_condition',0);
%%
br1_df=br_rvers(br1_df);
figure(443)
clf;
hold on
br1_df=br_contn(funcs_ss,br1_df,300);
br1_df=br_rvers(br1_df);
br1_df=br_contn(funcs_ss,br1_df,10000);
[nunst1_df,dom_df1,triv1_defect_df,br1_df.point]=GetStability(br1_df,'funcs',funcs_ss,...
    'exclude_trivial',true);%,'recompute',true);
chang_stb_df1=find(diff(nunst1_df));
% plotting one-par bif for symmetric POs in (df,max(u_A))-plan with stability
df1_xx=arrayfun(@(x)x.parameter(in.df),br1_df.point);
df1_pr_xx=arrayfun(@(x)x.parameter(in.PR),br1_df.point);
ymxx_df1=arrayfun(@(x)max(x.profile(1,:)),br1_df.point);
ymnn_df1=arrayfun(@(x)min(x.profile(1,:)),br1_df.point);
%% plot stability
figure(443)
clf;
plot(df1_xx(nunst1_df==0),ymxx_df1(nunst1_df==0),'ob',...
    df1_xx(nunst1_df>=1),ymxx_df1(nunst1_df>=1),'xk')
xlabel('df')
ylabel('max u_A')
grid on
%%
for i=1:length(br1_df.point)
    pp=br1_df.point(i);
    figure(66)
    plot(pp.mesh,pp.profile(1:2,:),'LineWidth',2)
end
%%  Branching off twoards non-symmetric solutions ( we vary df)
sbxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.5,6)');
poev1args={'usercond',{sbxsym},'initcond',{sbxsym}};
nspoev1args=addprefix('SetupPOEV1',poev1args);
[funcs_df1,nonsymper_df1,suc_v1]=SetupPsol(funcs_audi,br1_df,chang_stb_df1(2),'print_residual_info',1,...
   'outputfuncs',true,'branch_off','POEV1','contpar',in.df,...
    nspoev1args{:},'max_step',[in.PR,0.05; in.df,0.05; 0,0.01]);
figure(9993)
clf;
hold on
nonsymper_df1=br_contn(funcs_df1,nonsymper_df1,132);
%% Plloting stability for unsymmetric POs in (df,max(u_A))-plan
[nunst1_dfs,dom_df1s,triv1_defect_1dfs,nonsymper_df1.point]=GetStability(nonsymper_df1,'funcs',funcs_df1,...
    'exclude_trivial',true);%,'recompute',true);
chang_stb_df1=find(diff(nunst1_dfs));
df1_xxs=arrayfun(@(x)x.parameter(in.df),nonsymper_df1.point);
df1_xxs_pr=arrayfun(@(x)x.parameter(in.PR),nonsymper_df1.point);
ymxs_df1=arrayfun(@(x)max(x.profile(1,:)),nonsymper_df1.point);
ymns_df1=arrayfun(@(x)min(x.profile(1,:)),nonsymper_df1.point);
%% plot stability 
figure(4003)
clf;
hold on
plot(df1_xxs(nunst1_dfs==0),ymxs_df1(nunst1_dfs==0),'og',...
    df1_xxs(nunst1_dfs>=1),ymxs_df1(nunst1_dfs>=1),'xk')
plot(df1_xx(nunst1_df==0),ymxx_df1(nunst1_df==0),'ob',...
    df1_xx(nunst1_df>=1),ymxx_df1(nunst1_df>=1),'xk')
xlabel('d_f')
ylabel('max u_A')
title('varying d_f at fixed PR=25')
grid on
%
figure(4553)
clf;
hold on
plot(df1_xxs(nunst1_dfs==0),ymxs_df1(nunst1_dfs==0),'og',...
    df1_xxs(nunst1_dfs>=1),ymxs_df1(nunst1_dfs>=1),'xk')
br1_df=br_remove_extracolumns(br1_df);
nonsymper_df1=br_remove_extracolumns(nonsymper_df1);
%%
%save('unsymmetric_br_in_df_case_increased_D.mat')
%% Now from the above computation we can pick a non-symmetric PO with touching the threshold
% Then we continue in three parameters (rp,df,t) and plot the results in
% (rp,df)-space
c_A=[1,0,0,0,0,0];
smaxval=0.5;
uA_extrema1=arrayfun(@(p)dde_coll_roots(p,c_A,'diff',1)',nonsymper_df1.point,'uniformoutput',false);
ua_eval=@(p,t)c_A*dde_coll_eva(p.profile,p.mesh,t(:)',p.degree); % evaluate u_A at t in point p
sympomax1_ua=cellfun(@(p,t)max2(ua_eval(p,t)),num2cell(nonsymper_df1.point),uA_extrema1);
figure(9993)
clf;
plot(sympomax1_ua,'o')
grid on
%
[~,it_cross1]=min(abs(sympomax1_ua-(smaxval+1e-4)));
second_ua_peak1=3;
%
Rsym=[0,1,0,0,0,0;1,0,0,0,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
xdim=length(x0);
ine=in;
ine.t0=length(fieldnames(in))+1;
ine.val=ine.t0+1;
unsympo1=setfield(nonsymper_df1,'point',nonsymper_df1.point(it_cross1));
unsympo1.point.parameter([ine.t0,ine.val])=[uA_extrema1{it_cross1}(second_ua_peak1),smaxval];
max_cond=@(p,pref)dde_extreme_cond(p,c_A,ine.val,ine.t0);
%%
[mfuncs_df1,mbranch_df1,suc_max1]=ChangeBranchParameters(funcs_df1,unsympo1,1,...
    'contpar',[ine.PR,ine.df,ine.t0],...
    'usercond',{max_cond},'outputfuncs',true,...
    'print_residual_info',1,'max_step',[in.PR,0.1; in.df,0.05; 0,0.1]);
mbranch_df1.parameter.max_bound(3)=30;
figure(113)
%clf;
hold on
mbranch_df1=br_contn(mfuncs_df1,mbranch_df1,185);
mbranch_df1.parameter.max_step(4)=0.01;
mbranch_df1.parameter.max_step(5)=0.01;
mbranch_df1.parameter.max_step(6)=0.01;
%%
mbranch_df1.parameter.max_step(4)=0.08;
mbranch_df1.parameter.max_step(5)=0.05;
mbranch_df1.parameter.max_step(6)=0.08;
mbranch_df1=br_contn(mfuncs_df1,mbranch_df1,200);
%%
mbranch_df1=br_rvers(mbranch_df1);
mbranch_df1=br_contn(mfuncs_df1,mbranch_df1,200);
%% Compute and plot the stability
[mbranch_df1_wbifs,df1_tests,mm1_bifs,df1_bifind]=MonitorChange(mfuncs_df1,mbranch_df1,'print_residual_info',0);
nunst1_mm=GetStability(mbranch_df1_wbifs,'funcs',mfuncs_df1,...
    'exclude_trivial',true);
df1_m=arrayfun(@(x)x.parameter(in.df),mbranch_df1_wbifs.point);
pr1_m=arrayfun(@(x)x.parameter(in.PR),mbranch_df1_wbifs.point);
figure(33)
hold on %clf;
plot(pr1_m(nunst1_mm==0),df1_m(nunst1_mm==0),'g.',...
    pr1_m(nunst1_mm>=1),df1_m(nunst1_mm>=1),'y.','LineWidth',2)
grid on
%%
mbranch_df1_wbifs=br_remove_extracolumns(mbranch_df1_wbifs);
mbranch_df1=br_remove_extracolumns(mbranch_df1);
%
save('touching_theta_case_increased_TD_try2.mat')

%%












