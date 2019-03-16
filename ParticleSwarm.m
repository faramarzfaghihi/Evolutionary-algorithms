% Application of Particle swarm intelligence in minimization of functions
% writen by Faramarz Faghihi

clear all

N = 20;xl = -5;X_u = 5;CC2 = 2;CC1 = 2;
Theta_Max = 0.9;
Theta_Min  = 0.4;
i_Max = 200;
theta = Theta_Max;

V_old = zeros(2,N);
P_best = -100*ones(1,N);
G_best = -100;
count = 0; 

figure(1);clf;hold on
x=[-6:0.01:6];
y=[-6:0.01:6];
%--------------------------------------------------------------------
for i=1:N
    X(:,i) = xl+(X_u-xl)*rand(2,1); % generate random numbers 
    f(i) = X(1,i)^2+X(2,i)^2-(cos(2*pi*X(1,i))+cos(2*pi*X(2,i)));
end

F = 1./f;
for i=1:length(x)
    for j=1:length(y)
        f(j,i) = x(i)^2+y(j)^2-(cos(2*pi*x(i))+cos(2*pi*y(j)));
    end
end
contour(x,y,f,10)
colormap(gray)
xlabel('x_1')
ylabel('x_2')

for i=1:i_Max

    for j=1:N
        if F(j)>P_best(j)
            P_best(j) = F(j);
        end
    end
    
    if max(P_best)>G_best
        G_best = max(P_best);
        Temp = find(P_best==G_best,1);
        X_best = X(:,Temp);
        count = count+1;
    end
    
    V_New = theta*V_old+CC1*rand(2,N).*([P_best;P_best]-X)+CC2*rand(2,N).*(G_best*ones(2,N)-X);
    X = X+V_New;
    V_old = V_New;
    theta = Theta_Max-(Theta_Max-Theta_Min )*i/i_Max;
    Gbest_plot(i) = 1./G_best;
    
    for i=1:N
        f(i) = X(1,i)^2+X(2,i)^2-(cos(2*pi*X(1,i))+cos(2*pi*X(2,i)));
    end
    
    F = 1./f;
    plot(X_best(1,1),X_best(2,1),'k.')
    K = text(X_best(1,1),X_best(2,1)+.1,num2str(count));
    set(K,'FontSize',14)
    
end


figure(2)
plot(Gbest_plot,'k.')
xlabel('iteration number, i')
ylabel('G_{best}')

X_best
