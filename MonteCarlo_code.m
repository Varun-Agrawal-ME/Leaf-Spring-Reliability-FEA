%% Monte Carlo Simulation
function pf =MonteCarlo
L = [1149.2,91.1834232]; %mean and standard deviation of L respectively
b = [54.544,10.89660039]; %mean and standard deviation of b respectively
t = [7.209,1.762498227]; %mean and standard deviation of t respectively
W = [6000,500]; %mean and standard deviation of W respectively
E = 2.1e5; %Load applied at the center of the spring
n = 6; %total number of leaves
W = (3*W(1)*L(1))/(2*b(1)*t(1)^2*n); %initializing at the mean values
%Transformation from Lognormal to Normal for t
kesi_t = sqrt(log(1+(t(2)/t(1))^2));
lambda_t = log(t(1))-(kesi_t^2)/2;
nsamples = 140000; %number of samples
mu_Ln = 1149.2;
sigma_Ln = 91.1834232;
mu_bn = 54.544;
sigma_bn = 10.89660039;
mu_Wn = 6000;
sigma_Wn = 500;
R=rand(nsamples,1);
j=1;
for i=1:nsamples
    LS(i)=mu_Ln+sigma_Ln*norminv(R(i));
    bS(i)=mu_bn+sigma_bn*norminv(R(i));
    tS(i)=exp(lambda_t+kesi_t*norminv(R(i)));
    WS(i)=mu_Wn+sigma_Wn*norminv(R(i));
    G(i) = (3*WS(i).*LS(i))./(2.*bS(i).*tS(i).^2*n);
    if G(i)<125
        g(j)=G(i);
        j=j+1;
    else
    end
end
length(g)
pf=length(g)/nsamples; %probability of failure
Beta = abs(norminv(pf));
fprintf('The Reliability Index is equal to %s \n', Beta)
fprintf('The probability of failure is equal to %s', pf)
end