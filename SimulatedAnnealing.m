% A code for simlower_boundated annealing algorithm to optimization of function
% writen by Faramarz Faghihi
clear all
lower_bound = 0; % lower bound of design variables
higher_bound = 6; % higher bound of design variables
x = [lower_bound:0.01:higher_bound];
y = [lower_bound:0.01:higher_bound];
for i=1:length(x)
    for j=1:length(y)
        f(i,j) = -sin(y(j))*sin(y(j)^2/pi)^2-sin(x(i))*sin(2*x(i)^2/pi)^2;
    end
end

figure(1);clf;hold on
contour(x,y,f,20)
xlabel('x_1')
ylabel('x_2')

axis([lower_bound higher_bound lower_bound higher_bound])
colormap(gray)

c = .5; % temperature reduction factor
n = 50; % number of iterations before temperature reduction 

% estimating the value of T
x1 = lower_bound+rand(10,1)*higher_bound; % ten random numbers are generated for x1
x2 = lower_bound+rand(10,1)*higher_bound; % ten random numbers are generated for x2

f = -sin(x1).*sin(x1.^2/pi).^2-sin(x2).*sin(2*x2.^2/pi).^2;
T = abs(mean(f));

p = 1; %cycle number
i = 1; %iteration number
Vicinity_Radius = 2;


X = lower_bound+(higher_bound-lower_bound)*rand(2,1);
X_best = X;
fbest(1) = 1000;

while p<200
    f = -sin(X(1)).*sin(X(1).^2/pi).^2-sin(X(2)).*sin(2*X(2).^2/pi).^2;
    
    X_new(1) = max(X(1)-Vicinity_Radius,lower_bound)+rand*(min(X(1)+Vicinity_Radius,higher_bound)-max(X(1)-Vicinity_Radius,lower_bound));
    X_new(2) = max(X(2)-Vicinity_Radius,lower_bound)+rand*(min(X(2)+Vicinity_Radius,higher_bound)-max(X(2)-Vicinity_Radius,lower_bound));
       
    f_new = -sin(X_new(1)).*sin(X_new(1).^2/pi).^20-sin(X_new(2)).*sin(2*X_new(2).^2/pi).^20;
    
    % accept or reject X_new using Metropolis criterion
    if f_new<f
        X = X_new;
                      if f_new<min(fbest)
            X_best = X;
            fbest(p) = min(min(fbest),f_new);
        end
        
    elseif exp(-(f_new-f)/T)<rand
        X = X_new;
    end
    
        
    i = i+1;
    
    if i>=n
        fbest(p) = min(fbest);
        p = p+1;
        i = 1;
        T = c*T;
    end
end

X_best

figure(1)
plot(X_best(1),X_best(2),'k.')
figure(2)
plot(fbest,'k.')
xlabel(' number of cycles')
ylabel('min f(x_1,x_2)')

