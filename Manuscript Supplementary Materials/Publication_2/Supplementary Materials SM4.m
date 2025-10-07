fn = 'SM2 - Dataset series A1.xlsx';
global initial_conc data
sheetNames = sheetnames(fn);
% 37 is the most number of datapoint in a dataset
time = zeros(37,9);
data = zeros(37,9);
time(time==0)=NaN; 
data (data==0)=NaN;
initial_conc = [];

for k = 1:numel(sheetNames)
    [num,txt,raw] = xlsread(fn,sheetNames{k});
    initial_conc = [initial_conc; [num(4,2) num(5,2) num(6,2)] ];
    t = [];
    d = [];
    for j = 1:numel(num(9:end,1))
        time(j,k) = num(8+j,1);
        data(j,k) = num(8+j,2);
    end
    
end
%%
% K param: 1 := kcat , 2 := KIA , 3 := KMB, 4  := KMA, 5 := KIQ
% x concen: 1 := NAD+, 2 := HCO2 , 3 := NADH
%Formate dehydrogenase reaction : HCO2 + NAD+ <=> CO2 + NADH + H+
%Model only accounts for parameter fitting of HCO2, NAD+, and NADH
tic
K_fin = [0;0;0;0;0;0];

SSRs = [0]
for i = 1:1
    p0 = rand(6,1) * 1000; 
    [K_est, SSR,res,extflg, output] = lsqnonlin(@(p) errorFun(p,time,data), p0, zeros(1,6),[10000,10000,10000,10000,10000,10000] );
    tmp_SSR = SSR
    if SSR < 0.804*1.1
        K_fin = [K_fin K_est];
        SSRs = [SSRs SSR];
    end
end
toc
K_fin(:,1) = [];
mean(K_fin,2)
std(K_fin,[],2)
std(K_fin,[],2)./mean(K_fin,2)
save('results_10perc.mat')

function totalerr = errorFun(K, t, data)
    global initial_conc

    totalerr = [];
    for i = 1:9
       tt = t(:,i);
       tt(isnan(tt)) = []; %removes NaN
       dd = data(:,i);
       dd(isnan(dd)) = [];

       y_est = MMkinetics(K,tt,i);
       sd = dd./10;
       sd(1) = 1;
       err = (y_est - dd).^2;
       totalerr = [totalerr; [err]];
    end

    function [S] = MMkinetics(K,t,i)
    x0 = initial_conc(i,:);
    [T,Sv] = ode45(@DiffEq,t,x0);

        function dS = DiffEq(t,x)
            E = 0.022;
            xdot = zeros(3,1);
            xdot(1) = - K(1)*x(1)*x(2)/(K(2)*K(3) + K(3)*x(1) + K(4)*x(2) + x(1)*x(2) + K(2)*K(3)/K(5)*x(3) + K(4)/K(5)*x(2)*x(3) ) * E;
            xdot(2) = - K(1)*x(1)*x(2)/(K(2)*K(3) + K(3)*x(1) + K(4)*x(2) + x(1)*x(2) + K(2)*K(3)/K(5)*x(3) + K(4)/K(5)*x(2)*x(3) ) * E;
            xdot(3) = K(1)*x(1)*x(2)/(K(2)*K(3) + K(3)*x(1) + K(4)*x(2) + x(1)*x(2) + K(2)*K(3)/K(5)*x(3) + K(4)/K(5)*x(2)*x(3) ) * E - K(6)*x(3);
            dS = xdot;
        end
    S = Sv(:,3);
    end
end
