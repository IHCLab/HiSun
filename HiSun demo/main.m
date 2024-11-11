clear all; close all; clc;
%% Setting
load ALLINALL;
load wavelength_1997;
load Cuprite_1997_150_150_3
p= 9;
band1= [5:106]; band2= [115:151]; band3= [171:214]; banduse= [band1 band2 band3];
Y= X/10000; X=[];
[l,n]= size(Y);
[~,~,O_index]= RASF(Y,p,10/225);
Y(:,O_index)= [];
%% Algorithm
[M_est,X_est,time]= HiSun(Y,p);
%% Display
disp(['Time: ',num2str(time/60),' (min.)'])
sn_est= zeros(p,22500); sn_est(:,setdiff( [1:22500], O_index ))= X_est;
[~,A_scale]= rd_compare(M_est,ALLINALL,banduse);
A_scale= A_scale(:,p:-1:1);
range1= 1:length(band1);
range2= (range1(end)+1):(range1(end)+length(band2));
range3= (range2(end)+1):(range2(end)+length(band3));
step=0; figure;
for i=1:p
    L_min= min(A_scale(:,i))-0.1;
    L_max= max(A_scale(:,i))-L_min;
    A_test= A_scale(:,i)-ones(l,1)*L_min+step;
    step= step+L_max+0.1;
    plot(wavelength(range1),A_test(range1),wavelength(range2),A_test(range2),wavelength(range3),A_test(range3)); hold on
end
axis([0.4,2.5,0,step+0.4]);
sn_est_scale= sn_est./(max(sn_est')'*ones(1,n));
for i=1:p
    imnow= reshape(sn_est_scale(i,:),sqrt(n),sqrt(n));
    figure;
    imshow(imnow);title(i);
end