% A code for Artificial Bees colony optimization algorithms 
% writen by Faramarz Faghihi

clear all

lower_bound = 0;  % lower bound of design variables
higher_bound = 4; % higher bound of design variables

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

figure(2);clf;hold on
contour(x,y,f,20)
xlabel('x_1')
ylabel('x_2')

axis([lower_bound higher_bound lower_bound higher_bound])
colormap(gray)

SN = 30; % number of bees in colony
n = 2;   % number of varibles of the optimization problem
a = .1;

% initialization
x = lower_bound+(higher_bound-lower_bound)*rand(n,SN);

figure(1)
plot(x(1,:),x(2,:),'k.')

lim = 10;

for k=1:200
    
    % setting the position of employed bees
    for z=1:SN
        Temp1 = ceil(SN*rand);
        v(1,z) = x(1,z)+(-a+2*a*rand).*(x(1,z)-x(1,ceil(SN*rand)));
        while Temp1==z
            Temp1 = ceil(SN*rand);
            v(1,z) = x(1,z)+(-a+2*a*rand).*(x(1,z)-x(1,ceil(SN*rand)));
        end
        
        Temp1 = ceil(SN*rand);
        v(2,z) = x(2,z)+(-a+2*a*rand).*(x(2,z)-x(2,ceil(SN*rand)));
        while Temp1==z
            Temp1 = ceil(SN*rand);
            v(2,z) = x(2,z)+(-a+2*a*rand).*(x(2,z)-x(2,ceil(SN*rand)));
        end
    end
    
    % calclower_boundation of the fitness function
    Temp11 = -sin(v(1,:)).*sin(v(1,:).^2/pi).^2-sin(v(2,:)).*sin(2*v(2,:).^2/pi).^2;
   
    % onlooker bees
    for z=1:length(Temp11)
        if Temp11(z)>=0
            fitness(z) = 1/(1+Temp11(z));
        else
            fitness(z) = 1+abs(Temp11(z));
        end
    end
    
    pm = fitness/sum(fitness);
    
    % saving the best reslower_boundts obtained so far
    fbest(k) = min(-sin(x(1,:)).*sin(x(1,:).^2/pi).^2-sin(x(2,:)).*sin(2*x(2,:).^2/pi).^2);
    
    % calclower_boundation of the best solution obtained so far
    if fbest(k)<=min(fbest)
        x_best = v(:,find(fitness==max(fitness),1));
    else
        fbest(k) = min(fbest);
    end
    
    % implementation of the rolower_boundette-wheel
    for z = 1:SN 
        s = 0;
        Temp1 = rand;
        count = 0;
        while s<Temp1
            s = s+pm(count+1);
            count = count+1;
        end
        x(:,z) = v(:,count); %assignment of the new positions
    end
    
    % check if the bees find better solutions or not. If not, they sholower_boundd be
    % changed to scout bees
    if k>lim && all(diff(fbest(end-lim:end)))==0
        x = lower_bound+(higher_bound-lower_bound)*rand(n,SN);
    end
end

figure(2)
plot(x_best(1),x_best(2),'k.')

figure(3)
plot(fbest,'k')
xlabel('iteration number')
ylabel('min f(x_1,x_2)')

x_best
