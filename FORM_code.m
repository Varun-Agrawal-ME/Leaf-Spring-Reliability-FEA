%% FORM Method
%%%% Given limit state function is G() = (3*W*L)/(2*b*t^2*n)
% L is normally distributed, b is normally distributed,
% t is log-normally distributed, W is normally distributed
clc; clear; close all; format long;
%Given values
L = [1149.2,91.1834232]; %mean and standard deviation of L respectively
b = [54.544,10.89660039]; %mean and standard deviation of b respectively
t = [7.209,1.762498227]; %mean and standard deviation of t respectively
W = [6000,500]; %mean and standard deviation of W respectively
E = 2.1e5; %Load applied at the center of the spring
n = 6; %total number of leaves
G = (3*W(1)*L(1))/(2*b(1)*t(1)^2*n); %initializing at the mean values
%Transformation from Lognormal to Normal for t
kesi_t = sqrt(log(1+(t(2)/t(1))^2));
lambda_t = log(t(1))-0.5*(kesi_t^2);

count = 1; Betai = 1;
while G>125
    count=count +1
        %transformation to equivalent normal
        mu_Ln = 1149.2;
        sigma_Ln = 91.1834232;
        mu_bn = 54.544;
        sigma_bn = 10.89660039;
        mu_tn = t(1)*(1-log(t(1))+lambda_t);
        sigma_tn = t(1)*kesi_t;
        mu_Wn = 6000;
        sigma_Wn = 500;
        %transformation to reduced space
        L_pr = (L(1)-mu_Ln)/sigma_Ln;
        b_pr = (b(1)-mu_bn)/sigma_bn;
        t_pr = (t(1)-mu_tn)/sigma_tn;
        W_pr = (W(1)-mu_Wn)/sigma_Wn;
        %calculation of derivatives
        dg_L = 3*W(1)/(2*b(1)*t(1)^2*n);
        dg_b = -3*W(1)*L(1)/(2*b(1)^2*t(1)^2*n);
        dg_t = -3*W(1)*L(1)/(b(1)*t(1)^3*n);
        dg_W = 3*L(1)/(2*b(1)*t(1)^2*n);
        dgprMM = [dg_L*sigma_Ln;dg_b*sigma_bn;dg_t*sigma_tn;dg_W*sigma_Wn];
        %Gradient vector
        DG2 = dgprMM(1)^2+dgprMM(2)^2+dgprMM(3)^2+dgprMM(4)^2;
        %Calculate the new reduced coordinates
        L_prMn = (dgprMM(1)*L_pr + dgprMM(2)*b_pr + dgprMM(3)*t_pr + dgprMM(4)*W_pr - G)*dgprMM(1)/DG2;
        b_prMn = (dgprMM(1)*L_pr + dgprMM(2)*b_pr + dgprMM(3)*t_pr + dgprMM(4)*W_pr - G)*dgprMM(2)/DG2;
        t_prMn = (dgprMM(1)*L_pr + dgprMM(2)*b_pr + dgprMM(3)*t_pr + dgprMM(4)*W_pr - G)*dgprMM(3)/DG2;
        W_prMn = (dgprMM(1)*L_pr + dgprMM(2)*b_pr + dgprMM(3)*t_pr + dgprMM(4)*W_pr - G)*dgprMM(4)/DG2;
        %Transformation to original space
        L_new = L_prMn*sigma_Ln+mu_Ln;
        b_new = b_prMn*sigma_bn+mu_bn;
        t_new = t_prMn*sigma_tn+mu_tn;
        W_new = W_prMn*sigma_Wn+mu_Wn;
        %New value of limit state function
        G = (3*W_new*L_new)/(2*b_new*t_new^2*n)
        %Calculate Beta
        Beta = sqrt(L_prMn^2+b_prMn^2+t_prMn^2+W_prMn^2);
        Betadif = Beta - Betai;
        Betai = Beta;
        L(1) = L_new;
        b(1) = b_new;
        t(1) = t_new;
        W(1) = W_new;
        pf=1-normcdf(Betai);
end
fprintf('The final value of limit state function is equal to %s \n', G)
fprintf('The Reliability Index is equal to %s \n', Betai)
fprintf('The probability of failure is equal to %s', pf)