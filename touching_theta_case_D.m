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
[~,it]=min(abs(rp3_x-12));
[br3_df,sucp]=ChangeBranchParameters(funcs_ss,br_symmetry_wbifs(3),it,...
      'degree',6,'intervals',120,'contpar',in.df,'step',1e-3,'correc',true,parbds{:},...
      'extra_condition',true,'phase_condition',0);
%%
br3_df=br_rvers(br3_df);
figure(443)
clf;
hold on
br3_df=br_contn(funcs_ss,br3_df,300);
br3_df=br_rvers(br3_df);
br3_df=br_contn(funcs_ss,br3_df,10000);
[nunst3_df,dom_df3,triv3_defect_df,br3_df.point]=GetStability(br3_df,'funcs',funcs_ss,...
    'exclude_trivial',true);%,'recompute',true);
chang_stb_df3=find(diff(nunst3_df));
% plotting one-par bif for symmetric POs in (df,max(u_A))-plan with stability
df3_xx=arrayfun(@(x)x.parameter(in.df),br3_df.point);
df3_pr_xx=arrayfun(@(x)x.parameter(in.PR),br3_df.point);
ymxx_df3=arrayfun(@(x)max(x.profile(1,:)),br3_df.point);
ymnn_df3=arrayfun(@(x)min(x.profile(1,:)),br3_df.point);
%% plot stability
figure(443)
clf;
plot(df3_xx(nunst3_df==0),ymxx_df3(nunst3_df==0),'ob',...
    df3_xx(nunst3_df>=1),ymxx_df3(nunst3_df>=1),'xk')
xlabel('df')
ylabel('max u_A')
grid on
%%  Branching off twoards non-symmetric solutions ( we vary df)
sbxsym=@(p,pref)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0.1,0.5,6)');
poev1args={'usercond',{sbxsym},'initcond',{sbxsym}};
nspoev1args=addprefix('SetupPOEV1',poev1args);
[funcs_df3,nonsymper_df3,suc_v3]=SetupPsol(funcs_audi,br3_df,chang_stb_df3(2),'print_residual_info',1,...
   'outputfuncs',true,'branch_off','POEV1','contpar',in.df,...
    nspoev1args{:},'max_step',[in.PR,0.05; in.df,0.005; 0,0.001]);
figure(9993)
clf;
hold on
nonsymper_df3=br_contn(funcs_df3,nonsymper_df3,730);
%%
% nonsymper_df3=br_remove_extracolumns(nonsymper_df3);
% for i=1:length(nonsymper_df3.point)
% xp3=nonsymper_df3.point(i);
% figure(1000)
% clf;
% plot(xp3.mesh*xp3.period,xp3.profile(1:2,:),'LineWidth',2)
% grid on
% drawnow
% end
%% Plloting stability for unsymmetric POs in (df,max(u_A))-plan
[nunst3_dfs,dom_df3s,triv3_defect_3dfs,nonsymper_df3.point]=GetStability(nonsymper_df3,'funcs',funcs_df3,...
    'exclude_trivial',true);%,'recompute',true);
chang_stb_df3=find(diff(nunst3_dfs));
df3_xxs=arrayfun(@(x)x.parameter(in.df),nonsymper_df3.point);
df3_xxs_pr=arrayfun(@(x)x.parameter(in.PR),nonsymper_df3.point);
ymxs_df3=arrayfun(@(x)max(x.profile(1,:)),nonsymper_df3.point);
ymns_df3=arrayfun(@(x)min(x.profile(1,:)),nonsymper_df3.point);
%% plot stability 
figure(4003)
clf;
hold on
plot(df3_xxs(nunst3_dfs==0),ymxs_df3(nunst3_dfs==0),'og',...
    df3_xxs(nunst3_dfs>=1),ymxs_df3(nunst3_dfs>=1),'xk')
plot(df3_xx(nunst3_df==0),ymxx_df3(nunst3_df==0),'ob',...
    df3_xx(nunst3_df>=1),ymxx_df3(nunst3_df>=1),'xk')
xlabel('d_f')
ylabel('max u_A')
title('varying d_f at fixed PR=25')
grid on
%
figure(4553)
clf;
hold on
plot(df3_xxs(nunst3_dfs==0),ymxs_df3(nunst3_dfs==0),'og',...
    df3_xxs(nunst3_dfs>=1),ymxs_df3(nunst3_dfs>=1),'xk')
br3_df=br_remove_extracolumns(br3_df);
nonsymper_df3=br_remove_extracolumns(nonsymper_df3);
%%
%save('unsymmetric_br_in_df_case_increased_D.mat')
%% Now from the above computation we can pick a non-symmetric PO with touching the threshold
% Then we continue in three parameters (rp,df,t) and plot the results in
% (rp,df)-space
c_A=[1,0,0,0,0,0];
smaxval=0.5;
uA_extrema3=arrayfun(@(p)dde_coll_roots(p,c_A,'diff',1)',nonsymper_df3.point,'uniformoutput',false);
ua_eval=@(p,t)c_A*dde_coll_eva(p.profile,p.mesh,t(:)',p.degree); % evaluate u_A at t in point p
sympomax3_ua=cellfun(@(p,t)max2(ua_eval(p,t)),num2cell(nonsymper_df3.point),uA_extrema3);
figure(9993)
clf;
plot(sympomax3_ua,'o')
grid on
%
[~,it_cross3]=min(abs(sympomax3_ua-(smaxval+1e-4)));
second_ua_peak3=3;
%
Rsym=[0,1,0,0,0,0;1,0,0,0,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
xdim=length(x0);
ine=in;
ine.t0=length(fieldnames(in))+1;
ine.val=ine.t0+1;
unsympo3=setfield(nonsymper_df3,'point',nonsymper_df3.point(it_cross3));
unsympo3.point.parameter([ine.t0,ine.val])=[uA_extrema3{it_cross3}(second_ua_peak3),smaxval];
max_cond=@(p,pref)dde_extreme_cond(p,c_A,ine.val,ine.t0);
%%
[mfuncs_df3,mbranch_df3,suc_max3]=ChangeBranchParameters(funcs_df3,unsympo3,1,...
    'contpar',[ine.PR,ine.df,ine.t0],...
    'usercond',{max_cond},'outputfuncs',true,...
    'print_residual_info',1,'max_step',[in.PR,0.1; in.df,0.05; 0,0.1]);

mbranch_df3.parameter.max_bound(3)=30;

figure(113)
clf;
hold on
mbranch_df3=br_contn(mfuncs_df3,mbranch_df3,47);
mbranch_df3.parameter.max_step(4)=0.01;
mbranch_df3.parameter.max_step(5)=0.01;
mbranch_df3.parameter.max_step(6)=0.01;
mbranch_df3=br_contn(mfuncs_df3,mbranch_df3,90);
mbranch_df3.parameter.max_step(4)=0.05;
mbranch_df3.parameter.max_step(5)=0.05;
mbranch_df3.parameter.max_step(6)=0.05;
mbranch_df3=br_contn(mfuncs_df3,mbranch_df3,275);
%%
mbranch_df3=br_rvers(mbranch_df3);
mbranch_df3=br_contn(mfuncs_df3,mbranch_df3,200);
%% Compute and plot the stability
[mbranch_df3_wbifs,df3_tests,mm3_bifs,df3_bifind]=MonitorChange(mfuncs_df3,mbranch_df3,'print_residual_info',0);
nunst3_mm=GetStability(mbranch_df3_wbifs,'funcs',mfuncs_df3,...
    'exclude_trivial',true);
df3_m=arrayfun(@(x)x.parameter(in.df),mbranch_df3_wbifs.point);
pr3_m=arrayfun(@(x)x.parameter(in.PR),mbranch_df3_wbifs.point);
figure(33)
hold on %clf;
plot(pr3_m(nunst3_mm==0),df3_m(nunst3_mm==0),'g.',...
    pr3_m(nunst3_mm>=1),df3_m(nunst3_mm>=1),'y.','LineWidth',2)
grid on
%%
mbranch_df3_wbifs=br_remove_extracolumns(mbranch_df3_wbifs);
mbranch_df3=br_remove_extracolumns(mbranch_df3);
%
save('touching_theta_case_increased_D_try2.mat')

%%












