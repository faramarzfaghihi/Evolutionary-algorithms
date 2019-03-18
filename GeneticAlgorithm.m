% A simple code for Genetic algorithms for optimization of fucntions written by Faramarz Faghihi 

clear all

m = 50; % number of strings in mating pool
n = 2;  % number of variables
Q_bits = 11; % number of bits
l = 2*Q_bits;xu = 5;xl = -5;
Pc = 0.85;PR = 10;PC = ceil(Pc*m);
PM = m-(PR+PC);best(1) = 1000;
x = [-6:0.01:6]; 
y = [-6:0.01:6];

for i=1:length(x)
    for j=1:length(y)
        f(j,i) = (x(i)^2+y(j)-11)^2+(x(i)+y(j)^2-7)^2;
    end
end

figure(1);clf;hold on
contour(x,y,f,15)
plot(3,2,'xk')
plot(-2.8,3.13,'xk')
plot(-3.78,-3.28,'xk')
plot(3.58,-1.85,'xk')

xlabel('x_1')
ylabel('x_2')
colormap(gray)

clear x

% creating the initial mating pool
for k=1:m
    Mat_Pool(k,:) = round(rand(1,l));
end

for k=1:m
    x(k,1) = xl+(xu-xl)/(2^Q_bits-1)*bin2dec(num2str(Mat_Pool(k,1:Q_bits)));
    x(k,2) = xl+(xu-xl)/(2^Q_bits-1)*bin2dec(num2str(Mat_Pool(k,Q_bits+1:l)));
end

plot(x(:,1),x(:,2),'k.')

for v=1:100
    
    for k=1:m
        x(k,1) = xl+(xu-xl)/(2^Q_bits-1)*bin2dec(num2str(Mat_Pool(k,1:Q_bits)));
        x(k,2) = xl+(xu-xl)/(2^Q_bits-1)*bin2dec(num2str(Mat_Pool(k,Q_bits+1:l)));
        F(k,1) = 1/((x(k,1)^2+x(k,2)-11)^2+(x(k,1)+x(k,2)^2-7)^2);
    end

    
    % roulette wheel selection
    p = F/sum(F);
    P = cumsum(p);

    for k=1:PR
        pointr = find(P>rand,1);
        Mat_Pool_new(k,:) = Mat_Pool(pointr,:);
    end
    
     % implementation of the single point crossover
    for k=1:2:PC 
        parent1 = Mat_Pool(ceil(m*rand),:);
        parent2 = Mat_Pool(ceil(m*rand),:);
        pointr = ceil(l*rand);
        temp = parent1(pointr:l);
        parent1(pointr:l) = parent2(pointr:l);
        parent2(pointr:l) = temp;
        Mat_Pool_new(k+PR,:) = parent1;
        Mat_Pool_new(k+PR+1,:) = parent2;
    end

    % implementation of mutation
    for k=1:PM
        temp = Mat_Pool(ceil(m*rand),:);
        a = ceil(rand*l);
        if temp(a)==1
            temp(a) = 0;
        else
            temp(a) = 1;
        end
        Mat_Pool_new(PR+PC+k,:) = temp;
    end

    Mat_Pool = Mat_Pool_new;
    
    if 1/max(F)<min(best)
        best(v) = 1/max(F);
        best_location = find(F==max(F),1);
        X_opt(v,:) = x(best_location,:);
    else
        best(v) = min(best);
        X_opt(v,:) = X_opt(v-1,:);
    end
    
    
end

figure(2);clf;hold on

plot(x(:,1),x(:,2),'k.')

colormap(gray)

x = [-6:0.01:6];
y = [-6:0.01:6];
for i=1:length(x)
    for j=1:length(y)
        f(j,i) = (x(i)^2+y(j)-11)^2+(x(i)+y(j)^2-7)^2;
    end
end
contour(x,y,f,15)
plot(3,2,'xk')
plot(-2.8,3.13,'xk')
plot(-3.78,-3.28,'xk')
plot(3.58,-1.85,'xk')


xlabel('x_1')
ylabel('x_2')

figure(3);clf;hold on
plot(best,'k.')
xlabel('generation')
ylabel('min f(x_1,x_2)')

figure(4)
plot(X_opt(:,1),'k.')
xlabel('generation')
ylabel('optimum value of x_1')
figure(5)
plot(X_opt(:,2),'k.')
xlabel('generation')
ylabel('optimum value of x_2')

X_opt(end,:)
