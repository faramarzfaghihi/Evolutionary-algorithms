% Ant Colony Optimization to minimize functions  
% written by Faramarz Faghihi

clear all
X = [0:0.01:5];Y = [0:0.01:5];
n = 2;              % variables number
N = 50;             % the number of ants in the colony 
p = 1000;    
Rho = .5;Zeta = 20; % parameters values

%------------------------------------------------------------
for i=1:length(X)
    for j=1:length(Y)
        f(i,j) = -cos(Y(j))*sin(Y(j)^3/pi)^2-sin(X(i))*sin(2*X(i)^2/pi)^4;
        
    end
end
%------------------------------------------------------------

X1 = linspace(0,5,p);
X2 = linspace(0,5,p);
Tau = ones(2,length(X1));
Iteration_best_Fit(1) = 1000;
%--------------------------------------------------------
figure(1);clf;hold on
contour(X,Y,f,20)
xlabel('x1')
ylabel('x2')
axis([0 5 0 5])
colormap(gray)
%--------------------------------------------------------
figure(3);clf;hold on
contour(X,Y,f,20)
xlabel('x1')
ylabel('x2')
axis([0 5 0 5])
colormap(gray)
%--------------------------------------------------------

for k=1:100
        P = [Tau(1,1:end)/sum(Tau(1,1:end));Tau(2,1:end)/sum(Tau(2,1:end))]; %  probabilities matrix
        % implementing roulette-wheel
        cum_prob = [cumsum(P(1,1:end));cumsum(P(2,1:end))]; %cumulative probability
        cum_prob(:,1) = [0;0]; 
        for j=1:N
            Index = find(cum_prob(1,1:end)<=rand,1,'last');
            x1(1,j) = X1(Index);
            
            Index = find(cum_prob(2,1:end)<=rand,1,'last');
            x2(1,j) = X2(Index);
        end
        
        if k==1
            px = x1; py=x2;
        end
            
        f = -cos(Y(j))*sin(Y(j)^3/pi)^2-sin(X(i))*sin(2*X(i)^2/pi)^4;
        f_best = min(f);
        f_worst = max(f);
        Tau_old = (1-Rho)*Tau;
        Delta_Tau = abs(Zeta*f_best/(f_worst));
        
      
        if f_best<min(Iteration_best_Fit)
            Iteration_best_Fit(k) = f_best;
        else
            Iteration_best_Fit(k) = min(Iteration_best_Fit);
        end
            
                Index = find(f==f_best);
        
        for j=1:length(Index) % in general multiple ants use the best path
            i2 = find(X1==x1(Index(j)));
            j2 = find(X2==x2(Index(j)));
            Tau(1,i2) = Tau_old(1,i2)+Delta_Tau;
            Tau(2,j2) = Tau_old(2,j2)+Delta_Tau;
        end
               
end

x1_Best = X1(i2)
x2_Best = X2(j2)
figure(1)
plot(x1,x2,'k.')
%---------------------------------------------
figure(2)
plot(Iteration_best_Fit,'k.')
xlabel('Number of iterations')
ylabel('min f(x1,x2)')
%----------------------------------------------
figure(3)
plot(px,py,'k.')

