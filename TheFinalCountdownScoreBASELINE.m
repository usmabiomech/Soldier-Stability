function [slope]=TheFinalCountdownScoreBASELINE(RPE,avg_T_ang_std,D_avg_stances,D_std_stances,std_accel_1,std_accel_both,avg_peak_accel_1,avg_peak_accel_both)

% TheFinalCountdownScoreBASELINE(RPE,avg_T_ang_std,LyE_LF,D_LyE_TF,D_LyE_TS,D_LyE_SF,D_avg_stances,D_std_stances,std_accel_both,avg_peak_accel_both)

% make matrix for pc_base & slope w/ # variables
% slope=zeros(11,1);
slope=zeros(7,1);
% pc_base=zeros(11,1);
pc_base=zeros(7,1);

% Make time matrices
t1_matrix=1:length(avg_T_ang_std);
% t2_matrix=1:length(LyE_LF);
t4_matrix=1:length(D_avg_stances);
t6_matrix=1:length(avg_peak_accel_both);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS 1 - Trunk Angle Things
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[p_1,s_1]=polyfit(t1_matrix,avg_T_ang_std,1);
p_1_str=sprintf('y = %.4fx + %.3f',p_1(1),p_1(2));
R2_1=1-(s_1.normr/norm(avg_T_ang_std-mean(avg_T_ang_std)))^2;
R2_1_str=sprintf('R^2 = %.3f',R2_1);
f_1=polyval(p_1,t1_matrix);
pc_base(1)=(f_1(length(f_1))-f_1(1))/f_1(1);
slope(1)=RPE/pc_base(1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS 2 - Foot Acceleration Cycle
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % p_2=polyfit(t2_matrix,LyE_LF,1);
% mdl_2=fitlm(t2_matrix,LyE_LF);
% p_2=[mdl_2.Coefficients.Estimate(2),mdl_2.Coefficients.Estimate(1)];
% p_2_str=sprintf('y = %.4fx + %.3f',p_2(1),p_2(2));
% R2_2=mdl_2.Rsquared.Ordinary;
% R2_2_str=sprintf('R^2 = %.3f',R2_2);
% f_2=polyval(p_2,t2_matrix);
% pc_base(2)=(f_2(length(f_2))-f_2(1))/f_2(1);
% slope(2)=RPE/pc_base(2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS 3 - Change in LyE from T=>F, T=>S, and S=>F
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% % p_3_TF=polyfit(t2_matrix,D_LyE_TF,1);
% mdl_3_TF=fitlm(t2_matrix,D_LyE_TF);
% p_3_TF=[mdl_3_TF.Coefficients.Estimate(2),mdl_3_TF.Coefficients.Estimate(1)];
% p_3_TF_str=sprintf('y_{TF} = %.4fx + %.3f',p_3_TF(1),p_3_TF(2));
% R2_3_TF=mdl_3_TF.Rsquared.Ordinary;
% R2_3_TF_str=sprintf('R^2 = %.3f',R2_3_TF)
% f_3_TF=polyval(p_3_TF,t2_matrix);
% pc_base(3)=(f_3_TF(length(f_3_TF))-f_3_TF(1))/f_3_TF(1);
% slope(3)=RPE/pc_base(3);
% 
% % p_3_TS=polyfit(t2_matrix,D_LyE_TS,1);
% mdl_3_TS=fitlm(t2_matrix,D_LyE_TS);
% p_3_TS=[mdl_3_TS.Coefficients.Estimate(2),mdl_3_TS.Coefficients.Estimate(1)];
% p_3_TS_str=sprintf('y_{TS} = %.4fx + %.3f',p_3_TS(1),p_3_TS(2));
% R2_3_TS=mdl_3_TS.Rsquared.Ordinary;
% R2_3_TS_str=sprintf('R^2 = %.3f',R2_3_TS)
% f_3_TS=polyval(p_3_TS,t2_matrix);
% pc_base(4)=(f_3_TS(length(f_3_TS))-f_3_TS(1))/f_3_TS(1);
% slope(4)=RPE/pc_base(4);
% 
% % p_3_SF=polyfit(t2_matrix,D_LyE_SF,1);
% mdl_3_SF=fitlm(t2_matrix,D_LyE_SF);
% p_3_SF=[mdl_3_SF.Coefficients.Estimate(2),mdl_3_SF.Coefficients.Estimate(1)];
% p_3_SF_str=sprintf('y_{SF} = %.4fx + %.3f',p_3_SF(1),p_3_SF(2));
% R2_3_SF=mdl_3_SF.Rsquared.Ordinary;
% R2_3_SF_str=sprintf('R^2 = %.3f',R2_3_SF)
% f_3_SF=polyval(p_3_SF,t2_matrix);
% pc_base(5)=(f_3_SF(length(f_3_SF))-f_3_SF(1))/f_3_SF(1);
% slope(5)=RPE/pc_base(5);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESES 4 & 5 - Stance Times (AVG & SD)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[p_4,s_4]=polyfit(t4_matrix,D_avg_stances,1);
p_4_str=sprintf('y = %.4fx + %.3f',p_4(1),p_4(2));
R2_4=1-(s_4.normr/norm(D_avg_stances-mean(D_avg_stances)))^2;
R2_4_str=sprintf('R^2 = %.3f',R2_4);
f_4=polyval(p_4,t4_matrix);
% pc_base(6)=(f_4(length(f_4))-f_4(1))/f_4(1);
pc_base(2)=(f_4(length(f_4))-f_4(1))/f_4(1);
% slope(6)=RPE/pc_base(6);
slope(2)=RPE/pc_base(2);

[p_5,s_5]=polyfit(t4_matrix,D_std_stances,1);
p_5_str=sprintf('y = %.4fx + %.3f',p_5(1),p_5(2));
R2_5=1-(s_5.normr/norm(D_std_stances-mean(D_std_stances)))^2;
R2_5_str=sprintf('R^2 = %.3f',R2_5);
f_5=polyval(p_5,t4_matrix);
% pc_base(7)=(f_5(length(f_5))-f_5(1))/f_5(1);
% slope(7)=RPE/pc_base(7);
pc_base(3)=(f_5(length(f_5))-f_5(1))/f_5(1);
slope(3)=RPE/pc_base(3);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESES 6 & 7 - Total Foot Accelerations (AVG & SD)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[p_6,s_6]=polyfit(t6_matrix,avg_peak_accel_both,1);
p_6_str=sprintf('y_{Left+Right} = %.4fx + %.3f',p_6(1),p_6(2));
R2_6=1-(s_6.normr/norm(avg_peak_accel_both-mean(avg_peak_accel_both)))^2;
R2_6_str=sprintf('R^2 = %.3f',R2_6);
f_6=polyval(p_6,t6_matrix);
% pc_base(8)=(f_6(length(f_6))-f_6(1))/f_6(1);
% slope(8)=RPE/pc_base(8);
pc_base(4)=(f_6(length(f_6))-f_6(1))/f_6(1);
slope(4)=RPE/pc_base(4);

[p_6_1,s_6_1]=polyfit(t6_matrix,avg_peak_accel_1,1);
p_6_1_str=sprintf('y_{Left} = %.4fx + %.3f',p_6_1(1),p_6_1(2));
R2_6_1=1-(s_6_1.normr/norm(avg_peak_accel_1-mean(avg_peak_accel_1)))^2;
R2_6_1_str=sprintf('R^2 = %.3f',R2_6_1);
f_6_1=polyval(p_6_1,t6_matrix);
% pc_base(9)=(f_6_1(length(f_6_1))-f_6_1(1))/f_6_1(1);
% slope(9)=RPE/pc_base(9);
pc_base(5)=(f_6_1(length(f_6_1))-f_6_1(1))/f_6_1(1);
slope(5)=RPE/pc_base(5);

[p_7,s_7]=polyfit(t6_matrix, std_accel_both,1);
p_7_str=sprintf('y_{Left+Right} = %.4fx + %.3f',p_7(1),p_7(2));
R2_7=1-(s_7.normr/norm(std_accel_both-mean(std_accel_both)))^2;
R2_7_str=sprintf('R^2 = %.3f',R2_7);
f_7=polyval(p_7,t6_matrix);
% pc_base(10)=(f_7(length(f_7))-f_7(1))/f_7(1);
% slope(10)=RPE/pc_base(10);
pc_base(6)=(f_7(length(f_7))-f_7(1))/f_7(1);
slope(6)=RPE/pc_base(6);

[p_7_1,s_7_1]=polyfit(t6_matrix,std_accel_1,1);
p_7_1_str=sprintf('y_{Left} = %.4fx + %.3f',p_7_1(1),p_7_1(2));
R2_7_1=1-(s_7_1.normr/norm(std_accel_1-mean(std_accel_1)))^2;
R2_7_1_str=sprintf('R^2 = %.3f',R2_7_1);
f_7_1=polyval(p_7_1,t6_matrix);
% pc_base(11)=(f_7_1(length(f_7_1))-f_7_1(1))/f_7_1(1);
% slope(11)=RPE/pc_base(11);
pc_base(7)=(f_7_1(length(f_7_1))-f_7_1(1))/f_7_1(1);
slope(7)=RPE/pc_base(7);

end
