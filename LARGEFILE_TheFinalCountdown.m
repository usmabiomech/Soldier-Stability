function [avg_T_ang_std,D_avg_stances,D_std_stances,std_accel_1,std_accel_both,avg_peak_accel_1,avg_peak_accel_both] = TheFinalCountdown(filename)

% [avg_T_ang_std,LyE_LF,D_LyE_TF,D_LyE_TS,D_LyE_SF,D_avg_stances,D_std_stances,std_accel_both,avg_peak_accel_both]

% Takes in single file w/ info from all 4 dots => separates data into each
% dot (see address names below) => makes dots the same length => filters
% data => solves for (1) Trunk angle SD (2) Foot Acc LyE (3) Change b/t
% Foot/Shank/Trunk Acc LyE (4) Stance time AVG (5) Stance Time SD (6) Peak
% foot acc SD b/t L/R and (7) Peak foot acc AVG b/t L/R => writes data to
% (1) linearmetrics.csv and (2) nonlinearmetrics.csv (different time
% intervals therefore different sized matrices); currently outputs these individual
% metrics as well
% 
% file=uigetfile('*.csv');
% [timestamp, MAC, Eul_x, Eul_y, Eul_z, Acc_x, Acc_y, Acc_z] = chickenimport(filename);

% fid=fopen(filename); 
% datafile=textscan(fid, '%f %s %f %f %f %f %f %f','Delimiter',',','headerlines',6)
% fclose(fid);
 
% timestamp=datafile{1,1};
% MAC=string(datafile{1,2});
% Eul_x=datafile{1,3};
% Eul_y=datafile{1,4};
% Eul_z=datafile{1,5};
% Acc_x=datafile{1,6};
% Acc_y=datafile{1,7};
% Acc_z=datafile{1,8};

fid=fopen(filename);
datafile = fread(fid,'char');
fclose(fid)

data=cell(100000000,8);
newdata='';
m=1;
n=1;

for i=1:length(datafile)
    if datafile(i)==10
        m=m+1;
        n=1;
    elseif datafile(i)==44
        n=n+1;
    else 
        newdata=char(datafile(i));
        if i==1
            data{m,n}=newdata;
        else
            data{m,n}=horzcat(data{m,n},newdata);
        end
    end
end

% Make empty score array w/ number of score/pc_base/slope inputs
score=zeros(7,1);
pc_base=zeros(7,1);
slope=zeros(7,1);
T=zeros(10000000,6);
LS=zeros(10000000,6);
LF=zeros(10000000,6);
RF=zeros(10000000,6);
sig_x=zeros(1,1000000);
sig_y=zeros(1,1000000);
sig_z=zeros(1,1000000);
avg_L_stance=zeros(1,1000000);
avg_R_stance=zeros(1,1000000);
D_avg_stances=zeros(1,1000000);
std_L_stance=zeros(1,1000000);
std_R_stance=zeros(1,1000000);
D_std_stances=zeros(1,1000000);
avg_peak_accel_1=zeros(1,1000000);
avg_peak_accel_both=zeros(1,1000000);
std_accel_1=zeros(1,1000000);
std_accel_both=zeros(1,1000000);

% Find length of data, where first blank cell is 
i=7;
while i<=length(data)
    if ~isempty(data{i,1})
        i=i+1;
    else
        end_index=i-8; % -1 bc it is the index before; -7 bc starting @7
        break
    end
end

% Separate into variables; distinguish doubles and strings
timestamp=zeros(end_index,1);
MAC=cell(end_index,1);
Eul_x=zeros(end_index,1);
Eul_y=zeros(end_index,1);
Eul_z=zeros(end_index,1);
Acc_x=zeros(end_index,1);
Acc_y=zeros(end_index,1);
Acc_z=zeros(end_index,1);

for i=1:end_index
    timestamp(i)=str2double(horzcat(data{i+6,1}(2:end)));
    MAC{i}=data{i+6,2};
    Eul_x(i)=str2double(horzcat(data{i+6,3}));
    Eul_y(i)=str2double(horzcat(data{i+6,4}));
    Eul_z(i)=str2double(horzcat(data{i+6,5}));
    Acc_x(i)=str2double(horzcat(data{i+6,6}));
    Acc_y(i)=str2double(horzcat(data{i+6,7}));
    Acc_z(i)=str2double(horzcat(data{i+6,8}(1:end-2)));
end

% Dot 1 = Trunk = d4:22:cd:00:03:82
% Dot 2 = Left Shank = d4:22:cd:00:03:79
% Dot 3 = Left Foot = d4:22:cd:00:03:8a
% Dot 4 = Right Shank = d4:22:cd:00:03:76
% DOt 5 = Right Foot = d4:22:cd:00:03:92
n_t=1; n_ls=1; n_lf=1; n_rs=1; n_rf=1;
for i=1:length(timestamp)
    if MAC{i}=='d4:22:cd:00:03:82'
        T(n_t,:)=[Eul_x(i) Eul_y(i) Eul_z(i) Acc_x(i) Acc_y(i) Acc_z(i)];
        n_t=n_t+1;
    elseif MAC{i}=='d4:22:cd:00:03:79'
        LS(n_ls,:)=[Eul_x(i) Eul_y(i) Eul_z(i) Acc_x(i) Acc_y(i) Acc_z(i)];
        n_ls=n_ls+1;
    elseif MAC{i}=='d4:22:cd:00:03:8a'
        LF(n_lf,:)=[Eul_x(i) Eul_y(i) Eul_z(i) Acc_x(i) Acc_y(i) Acc_z(i)];
        n_lf=n_lf+1;
%     elseif data_txt{i,1}=='d4:22:cd:00:03:95'
%         RS(n_rs,:)=data_num(i,3:8);
%         n_rs=n_rs+1;
    elseif MAC{i}=='d4:22:cd:00:03:76'
        RF(n_rf,:)=[Eul_x(i) Eul_y(i) Eul_z(i) Acc_x(i) Acc_y(i) Acc_z(i)];
        n_rf=n_rf+1;
    end
end

% Cut off zeros from end of T, LS, LF, RF matrices

T=zeroKillerChicken(T);
LS=zeroKillerChicken(LS);
LF=zeroKillerChicken(LF);
RF=zeroKillerChicken(RF);

% Make all the same length
cap=min([length(T(:,1)) length(LS(:,1)) length(LF(:,1)) length(RF(:,1))]);
T=T(1:cap,:);
LS=LS(1:cap,:);
LF=LF(1:cap,:);
RF=RF(1:cap,:);
    
% Columns: (1)Euler_X (2)Euler_Y (3)Euler_Z (4)FreeAcc_x (5)FreeAcc_y
% (5)FreeAcc_Z

% Butterworth filtering
fc=10;%cutoff frequency
fs=60;%DOT sampled frequency
[b,a]=butter(4,fc/(fs/2));%coefficients for low-pass 4th order butterworth filter
%a cutoff frequency of 10
for i=1:length(T(1,:))
    T(:,i)=filter(b,a,T(:,i));
    LS(:,i)=filter(b,a,LS(:,i));
    LF(:,i)=filter(b,a,LF(:,i));
    RF(:,i)=filter(b,a,RF(:,i));
end
time=0:(length(T(:,1))-1); %[Hz]

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS 1 - Trunk Angle Things
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t1=1; % length of time interval [min]
t1_matrix=1:(floor(length(time)/fs/60/t1));
j1=1;
k1=fs*60*t1;
for n=1:(floor(length(time)/fs/60/t1)) %n=number of full n-minute intervals
    sig_x(n)=std(T(j1:k1,4));
    sig_y(n)=std(T(j1:k1,5));
    sig_z(n)=std(T(j1:k1,6));
    j1=k1+1;
    k1=j1+(fs*60*t1-1);
end
std_T=vertcat(sig_x,sig_y,sig_z);
avg_T_ang_std=mean(std_T);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS 2 - Foot Acceleration Cycle
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% t2=5; % length of time interval [min]
% t2_matrix=1:(floor(length(time)/fs/60/t2));
% j2=1;
% k2=fs*60*t2;
% 
% LF_acc_tot=sqrt(LF(:,4).^2+LF(:,5).^2+LF(:,6).^2);
% T_acc_tot=sqrt(T(:,4).^2+T(:,5).^2+T(:,6).^2);
% LS_acc_tot=sqrt(LS(:,4).^2+LS(:,5).^2+LS(:,6).^2);
% dim=3;
% eRange=400;
% 
% for n=1:(floor(length(time)/fs/60/t2)) %n=number of full n-minute intervals
%     [~,lag]=phaseSpaceReconstruction(LF_acc_tot(j2:k2),[],dim);
%     LyE_LF(n)=lyapunovExponent(LF_acc_tot(j2:k2),fs,lag,dim,'ExpansionRange',eRange);
%     [~,lag]=phaseSpaceReconstruction(T_acc_tot(j2:k2),[],dim);
%     LyE_T(n)=lyapunovExponent(T_acc_tot(j2:k2),fs,lag,dim,'ExpansionRange',eRange);
%     [~,lag]=phaseSpaceReconstruction(LS_acc_tot(j2:k2),[],dim);
%     LyE_LS(n)=lyapunovExponent(LS_acc_tot(j2:k2),fs,lag,dim,'ExpansionRange',eRange);
%     j2=k2+1;
%     k2=j2+(fs*60*t1-1);
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS 3 - Change in LyE from T=>F, T=>S, and S=>F
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% D_LyE_TF=abs(LyE_T-LyE_LF); % change in LyE from Trunk to Shank
% D_LyE_TS=abs(LyE_T-LyE_LS);
% D_LyE_SF=abs(LyE_LS-LyE_LF);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESES 4 & 5 - Stance Times (AVG & SD)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t4=1; % length of time interval [min]
t4_matrix=1:(floor(length(time)/fs/60/t4));
j4=1;
k4=fs*60*t4;
stancetimes_L = 0;
threshold = 6;
min_length = 15;
stancetimes_R = 0;
threshold_R = 6;
min_length_R = 15;
ii = 1;
for n=1:(floor(length(time)/fs/60/t4))
    accel_x_L = LF(j4:k4,4);
    accel_x_R = RF(j4:k4,4);
    i = 1;
    while i <= length(accel_x_L)
        if accel_x_L(i) < threshold
            c = 0;
            while accel_x_L(i+c) < threshold
                c = c+1;
                if i+c >length(accel_x_L)
                    break
                end
            end
            if c > min_length
                stancetimes_L = [stancetimes_L,c];
                i = i+c;
            end
        else
            i = i+1;
        end
        i = i+1;
    end
    
    while ii <= length(accel_x_R)
        if accel_x_R(ii) < threshold_R
            cc = 0;
            while accel_x_R(ii+cc) < threshold_R
                cc = cc+1;
                if ii+cc >length(accel_x_R)
                    break
                end
            end
            if cc > min_length_R
                stancetimes_R = [stancetimes_R,cc];
                ii = ii+cc;
            end
        else
            ii = ii+1;
        end
        ii = ii+1;
    end

    i = 1;
    while i <= length(stancetimes_L)
        if stancetimes_L(i) >200
            stancetimes_L=[stancetimes_L(1:i-1),stancetimes_L(i+1:end)];
        else
            i = i+1;
        end
    end

    i = 1;
    while i <= length(stancetimes_R)
        if stancetimes_R(i) >200
            stancetimes_R=[stancetimes_R(1:i-1),stancetimes_R(i+1:end)];
        else
            i = i+1;
        end
    end
    avg_L_stance(n) = mean(stancetimes_L);
    avg_R_stance(n) = mean(stancetimes_R);
    D_avg_stances(n) = abs(avg_L_stance(n) - avg_R_stance(n));

    std_L_stance(n) = std(stancetimes_L);
    std_R_stance(n) = std(stancetimes_R);
    D_std_stances(n) = abs(std_L_stance(n) - std_R_stance(n));
    
    j4=k4+1;
    k4=j4+(fs*60*t4-1);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESES 6 & 7 - Total Foot Accelerations (AVG & SD)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t6=1; % length of time interval [min]
t6_matrix=1:(floor(length(time)/fs/60/t6));
j6=1;
k6=fs*60*t6;
for n=1:(floor(length(time)/fs/60/t6))
    total_LF = sqrt(LF(j6:k6,4).^2 + LF(j6:k6,5).^2 + LF(j6:k6,6).^2) - 9.8;
    max_LF = islocalmax(total_LF,'MinSeparation',fs);
    loc_max_LF = find(max_LF ==1);
    total_RF = sqrt(RF(j6:k6,4).^2 + RF(j6:k6,5).^2 + RF(j6:k6,6).^2) - 9.8;
    max_RF = islocalmax(total_RF,'MinSeparation',fs);
    loc_max_RF = find(max_RF ==1);

    avg_L = mean(total_LF(loc_max_LF));
    err_L = std(total_LF(loc_max_LF));
    avg_R = mean(total_RF(loc_max_RF));
    err_R = std(total_RF(loc_max_RF));

    avg_both = abs(avg_L+avg_R);
    avg_peak_accel_1(n) = avg_L;
    avg_peak_accel_both(n) = avg_both;

    std_both = abs(err_L+err_R);
    std_accel_1(n) = err_L;
    std_accel_both(n) = std_both;
    
    j6=k6+1;
    k6=j6+(fs*60*t6-1);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Cut off zeros from end of sig_x, sig_y, sig_z, avg_L_stance, avg_R_stance,
% D_avg_stances, std_L_stance, std_R_stance, D_std_stances, avg_peak_accel_1, 
% avg_peak_accel_both, std_accel_1, std_accel_both

avg_T_ang_std=zeroKillerChicken2(avg_T_ang_std);
D_avg_stances=zeroKillerChicken2(D_avg_stances);
D_std_stances=zeroKillerChicken2(D_std_stances);
std_accel_1=zeroKillerChicken2(std_accel_1);
std_accel_both=zeroKillerChicken2(std_accel_both);
avg_peak_accel_1=zeroKillerChicken2(avg_peak_accel_1);
avg_peak_accel_both=zeroKillerChicken2(avg_peak_accel_both);

end