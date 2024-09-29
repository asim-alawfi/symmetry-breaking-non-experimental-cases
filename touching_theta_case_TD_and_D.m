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
[~,it]=min(abs(rp2_x-10.1));
[br2_df,sucp]=ChangeBranchParameters(funcs_ss,br_symmetry_wbifs(2),it,...
      'degree',6,'intervals',120,'contpar',in.df,'step',1e-3,'correc',true,parbds{:},...
      'extra_condition',true,'phase_condition',0);
%%
br2_df=br_rvers(br2_df);
figure(443)
clf;
hold on
br2_df=br_contn(funcs_ss,br2_df,300);
br2_df=br_rvers(br2_df);
br2_df=br_contn(funcs_ss,br2_df,10000);
[nunst2_df,dom_df2,triv2_defect_df,br2_df.point]=GetStability(br2_df,'funcs',funcs_ss,...
    'exclude_trivial',true);%,'recompute',true);
chang_stb_df2=find(diff(nunst2_df));
% plotting one-par bif for symmetric POs in (df,max(u_A))-plan with stability
df2_xx=arrayfun(@(x)x.parameter(in.df),br2_df.point);
df2_pr_xx=arrayfun(@(x)x.parameter(in.PR),br2_df.point);
ymxx_df2=arrayfun(@(x)max(x.profile(1,:)),br2_df.point);
ymnn_df2=arrayfun(@(x)min(x.profile(1,:)),br2_df.point);
%% plot stability
figure(443)
clf;
plot(df2_xx(nunst2_df==0),ymxx_df2(nunst2_df==0),'ob',...
    df2_xx(nunst2_df>=1),ymxx_df2(nunst2_df>=1),'xk')
xlabel('df')
ylabel('max u_A')
grid on
%%
for i=1:length(br2_df.point)
    pp=br2_df.point(i);
    figure(66)
    plot(pp.mesh,pp.profile(1:2,:),'LineWidth',2)
    yline(0.5,'k--')
end
%%  Branching off twoards non-symmetric solutions ( we vary df)
sbxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.5,6)');
poev1args={'usercond',{sbxsym},'initcond',{sbxsym}};
nspoev1args=addprefix('SetupPOEV1',poev1args);
[funcs_df2,nonsymper_df2,suc_v2]=SetupPsol(funcs_audi,br2_df,chang_stb_df2(2),'print_residual_info',1,...
   'outputfuncs',true,'branch_off','POEV1','contpar',in.df,...
    nspoev1args{:},'max_step',[in.PR,0.05; in.df,0.05; 0,0.01]);
figure(9993)
clf;
hold on
nonsymper_df2=br_contn(funcs_df2,nonsymper_df2,132);
%% Plloting stability for unsymmetric POs in (df,max(u_A))-plan
[nunst2_dfs,dom_df2s,triv2_defect_2dfs,nonsymper_df2.point]=GetStability(nonsymper_df2,'funcs',funcs_df2,...
    'exclude_trivial',true);%,'recompute',true);
chang_stb_df2=find(diff(nunst2_dfs));
df2_xxs=arrayfun(@(x)x.parameter(in.df),nonsymper_df2.point);
df2_xxs_pr=arrayfun(@(x)x.parameter(in.PR),nonsymper_df2.point);
ymxs_df2=arrayfun(@(x)max(x.profile(1,:)),nonsymper_df2.point);
ymns_df2=arrayfun(@(x)min(x.profile(1,:)),nonsymper_df2.point);
%% plot stability 
figure(4003)
clf;
hold on
plot(df2_xxs(nunst2_dfs==0),ymxs_df2(nunst2_dfs==0),'og',...
    df2_xxs(nunst2_dfs>=1),ymxs_df2(nunst2_dfs>=1),'xk')
plot(df2_xx(nunst2_df==0),ymxx_df2(nunst2_df==0),'ob',...
    df2_xx(nunst2_df>=1),ymxx_df2(nunst2_df>=1),'xk')
xlabel('d_f')
ylabel('max u_A')
title('varying d_f at fixed PR=25')
grid on
%
figure(4553)
clf;
hold on
plot(df2_xxs(nunst2_dfs==0),ymxs_df2(nunst2_dfs==0),'og',...
    df2_xxs(nunst2_dfs>=1),ymxs_df2(nunst2_dfs>=1),'xk')
br2_df=br_remove_extracolumns(br2_df);
nonsymper_df2=br_remove_extracolumns(nonsymper_df2);
%%
%save('unsymmetric_br_in_df_case_increased_D.mat')
%% Now from the above computation we can pick a non-symmetric PO with touching the threshold
% Then we continue in three parameters (rp,df,t) and plot the results in
% (rp,df)-space
c_A=[1,0,0,0,0,0];
smaxval=0.5;
uA_extrema2=arrayfun(@(p)dde_coll_roots(p,c_A,'diff',1)',nonsymper_df2.point,'uniformoutput',false);
ua_eval=@(p,t)c_A*dde_coll_eva(p.profile,p.mesh,t(:)',p.degree); % evaluate u_A at t in point p
sympomax2_ua=cellfun(@(p,t)max2(ua_eval(p,t)),num2cell(nonsymper_df2.point),uA_extrema2);
figure(9993)
clf;
plot(sympomax2_ua,'o')
grid on
%
[~,it_cross2]=min(abs(sympomax2_ua-(smaxval+1e-4)));
second_ua_peak2=3;
%
Rsym=[0,1,0,0,0,0;1,0,0,0,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
xdim=length(x0);
ine=in;
ine.t0=length(fieldnames(in))+1;
ine.val=ine.t0+1;
unsympo2=setfield(nonsymper_df2,'point',nonsymper_df2.point(it_cross2));
unsympo2.point.parameter([ine.t0,ine.val])=[uA_extrema2{it_cross2}(second_ua_peak2),smaxval];
max_cond=@(p,pref)dde_extreme_cond(p,c_A,ine.val,ine.t0);
%%
[mfuncs_df2,mbranch_df2,suc_max2]=ChangeBranchParameters(funcs_df2,unsympo2,1,...
    'contpar',[ine.PR,ine.df,ine.t0],...
    'usercond',{max_cond},'outputfuncs',true,...
    'print_residual_info',1,'max_step',[in.PR,0.08; in.df,0.05; 0,0.08]);
%
mbranch_df2.parameter.max_bound(3)=30;
figure(113)
clf;
hold on
mbranch_df2=br_contn(mfuncs_df2,mbranch_df2,70);
%%
mbranch_df2.parameter.max_step(4)=0.01;
mbranch_df2.parameter.max_step(5)=0.01;
mbranch_df2.parameter.max_step(6)=0.01;
mbranch_df2=br_contn(mfuncs_df2,mbranch_df2,200);
mbranch_df2.parameter.max_step(4)=0.08;
mbranch_df2.parameter.max_step(5)=0.05;
mbranch_df2.parameter.max_step(6)=0.08;
mbranch_df2=br_contn(mfuncs_df2,mbranch_df2,300);
%%
mbranch_df2=br_rvers(mbranch_df2);
mbranch_df2=br_contn(mfuncs_df2,mbranch_df2,200);
%% Compute and plot the stability
[mbranch_df2_wbifs,df2_tests,mm2_bifs,df2_bifind]=MonitorChange(mfuncs_df2,mbranch_df2,'print_residual_info',0);
nunst2_mm=GetStability(mbranch_df2_wbifs,'funcs',mfuncs_df2,...
    'exclude_trivial',true);
df2_m=arrayfun(@(x)x.parameter(in.df),mbranch_df2_wbifs.point);
pr2_m=arrayfun(@(x)x.parameter(in.PR),mbranch_df2_wbifs.point);
figure(33)
hold on %clf;
plot(pr2_m(nunst2_mm==0),df2_m(nunst2_mm==0),'g.',...
    pr2_m(nunst2_mm>=1),df2_m(nunst2_mm>=1),'y.','LineWidth',2)
grid on
%%
mbranch_df2_wbifs=br_remove_extracolumns(mbranch_df2_wbifs);
mbranch_df2=br_remove_extracolumns(mbranch_df2);
%
save('touching_theta_case_increased_TD_and_D_try2.mat')

%%












