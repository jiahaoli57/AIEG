% function [daily_return, total_return] = AIEG_run(data)
%{
This file is the run core for the Competitive online strategy based on 
improved exponential gradient expert and aggregating method (AIEG).

For any usage of this function, the following papers should be cited as
reference:

[1] Yong Zhang, Jiahao Li, Xingyu Yang, and Jianliang Zhang. "Competitive 
online strategy based on improved exponential gradient expert and 
aggregating method" Computational Economics, 2023.

Inputs:
data                      -data with price relative sequences

Outputs:
daily_return              -daily wealths
total_return              -total wealths
%}

%% Parameter Setting
tc = 0; % transaction cost rate
eta_min = 0.01;
step = 0.01;
eta_max = 0.2;
w = 5;

%% Variables Inital
[T,N] = size(data);
b = zeros(T,N);
daily_return = zeros(T,1);

s = cell(1,N);
for ns = 1:N
    s{1,ns} = 0;
end

num_eta=0;
for eta = eta_min:step:eta_max
    num_eta = num_eta+1;
end

S = cell(1,N);
for nS = 1:N
    S{1,nS} = ones(num_eta,1);
end

e = cell(T,num_eta);
h = ones(1,N)/N;
for i = 1:num_eta
    e{1,i} = h;
end

%% Calculate the close prices
for t = 1:T
    if t == 1
        data_close(t,:) = data(t,:);
    else
        data_close(t,:) = data(t,:).*data_close(t-1,:);
  end
end

%% Main
for t = 1:T
    if t==1
        b(1,:) = ones(1,N)/N;
        daily_return(t,1) = b(t,:)*data(t,:)';
        
%         daliy_exp_r(:,:,t) = b(:,:)*data(t,:)';
%         exp_cumres(:,:,t) = daliy_exp_r(:,:,t);
    else
        t1=t;
        if t1<w+2                                   
            x_t1(t1,:) = data(t1-1,:);
        else
            x_t1(t1,:) = l1median_VaZh_z(data_close((t1-w):(t1-1),:))./data_close(t1-1,:);
        end

        k=0;                                                             
        for eta = eta_min:step:eta_max
            k = k+1;
            ff = e{t-1,k};
            Z = ff.*exp(eta*x_t1(t1,:)/(ff*x_t1(t1,:)'))*ones(N,1);
            f = ff.*exp(eta*x_t1(t1,:)/(ff*x_t1(t1,:)'));
            f = f/Z; 
            e{t,k} = f;


            exp_h(t-1,:)=data(t-1,:).*e{t-1,k};
            exp_hat(t-1,:)=exp_h(t-1,:)/sum(exp_h(t-1,:));
            exp_diff(t,:)=sum(abs(e{t,k}-exp_hat(t-1,:)));  

            for n1 = 1:N
                 S{1,n1}(k,1) = S{1,n1}(k,1)*(ff*data(t-1,:)'*(1-(tc)*exp_diff(t-1,:)));
            end

            s{1,N} = s{1,N}+S{1,N}(k,1)^(1/sqrt(t));
            for n2 = 1:N-1
                s{1,n2} = s{1,n2}+S{1,n2}(k,1)^(1/sqrt(t))*f(:,n2);
            end
        end
        
        SUM_N = 0;
        for n3 = 1:N
            if n3 < N
                b(t,n3) = s{1,n3}/s{1,N};
                SUM_N = SUM_N+s{1,n3};
            else
                b(t,n3) = (s{1,N}-SUM_N)/s{1,N};
            end
        end

        b_h(t-1,:)=data(t-1,:).*b(t-1,:);
        b_hat(t-1,:)=b_h(t-1,:)/sum(b_h(t-1,:));
        diff(t,:)=sum(abs(b(t,:)-b_hat(t-1,:)));
     

        daily_return(t) = (data(t,:)*b(t,:)'*(1-(tc)*diff(t,:)));
    end
end
total_return = cumprod(daily_return);
