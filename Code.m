%% Real-Coded Genetic Algorithem 
clearvars; clc;
%% Determine constants 
Pc = 0.9;                          %crossover population
Pm = 0.05;                         %mutation rate 
gen = 300;                         %number of generation required
kmin = 1; kmax = 50;               %Gain limits 
T1min = 0.1; T1max = 0.5;          %Time constant 1 limits
T2min = 0.1; T2max = 0.1;          %Time constant 2 limits
limits = [kmin, T1min, T2max; kmax, T1max, T2max]; %limits mat. of par.
alpha = 1;                       %To implement BLX-alpha
beta = 1;                          %B-operator of non-uniform mutation 
npop = 50;                         %Number of population

%% RCGA for 10 times
for run=1:1:10
    % Generate initial population
    k = kmin + (kmax-kmin)*rand(npop,1);
    T1 = T1min + (T1max-T1min)*rand(npop,1);
    T2 = T2min + (T2max-T2min)*rand(npop,1);
    old_pop = [k,T1,T2];

    % Main loop of RCGA
    for iter=1:1:gen
        new_pop = zeros(npop,3);
        % Tournement Selection
        for i = 1:1:npop
            k = ceil(npop*rand(1,1)); %choose the first chromoshom location
            j = ceil(npop*rand(1,1)); %chosse the second chromoshom location
            a = ff(old_pop(k,:));
            b = ff(old_pop(j,:));
                if a >= b;            %choose the fittest chromoshom
                    new_pop(i,:) = old_pop(k,:); 
                else 
                    new_pop(i,:) = old_pop(j,:);
                end
        end
        % BLX-alpha Crossover
        temp_pop=new_pop;
        for i = 1:1:npop/2  
            cros_prop = rand;
            if cros_prop <= Pc
                BLX = [min(temp_pop(2*i-1,1),temp_pop(2*i,1))-abs(temp_pop(2*i-1,1)...
                    -temp_pop(2*i,1))*alpha , min(temp_pop(2*i-1,2),temp_pop(2*i,2))-...
                    abs(temp_pop(2*i-1,2)-temp_pop(2*i,2))*alpha , min(temp_pop(2*i-1,3),...
                    temp_pop(2*i,3))-abs(temp_pop(2*i-1,3)-temp_pop(2*i,3))*alpha; ...
                    max(temp_pop(2*i-1,1),temp_pop(2*i,1))+abs(temp_pop(2*i-1,1)...
                    -temp_pop(2*i,1))*alpha , max(temp_pop(2*i-1,2),temp_pop(2*i,2))+...
                    abs(temp_pop(2*i-1,2)-temp_pop(2*i,2))*alpha , max(temp_pop(2*i-1,3),...
                    temp_pop(2*i,3))+abs(temp_pop(2*i-1,3)-temp_pop(2*i,3))*alpha]; 
                for j = 1:1:3
                    new_pop(2*i-1,j) = BLX(1,j) + [BLX(2,j)-BLX(1,j)]*rand;
                    new_pop(2*i,j) = BLX(1,j) + [BLX(2,j)-BLX(1,j)]*rand;
                    if new_pop(2*i-1,j) < limits(1,j) || new_pop(2*i-1,j) > limits(2,j)
                        new_pop(2*i-1,j) = min(new_pop(2*i-1,j),limits(2,j));
                        new_pop(2*i-1,j) = max(new_pop(2*i-1,j),limits(1,j));
                    end
                    if new_pop(2*i,j) < limits(1,j) || new_pop(2*i,j) > limits(2,j)
                        new_pop(2*i,j) = min(new_pop(2*i,j),limits(2,j));
                        new_pop(2*i,j) = max(new_pop(2*i,j),limits(1,j));
                    end                           
                end       
            end
        end
        % Non-unifrom Mutation
        for i=1:1:npop
            for j = 1:1:3
                mu_prop = rand;
                if mu_prop <= Pm
                    delta = (1-rand(1,1)^((1-iter/gen)^beta));
                    tao = round(rand);
                    if tao == 1
                        new_pop(i,j) = new_pop(i,j) - delta*(new_pop(i,j)-limits(1,j));
                    else 
                        new_pop(i,j) = new_pop(i,j) + delta*(limits(2,j)-new_pop(i,j));
                    end
                end
                    if new_pop(i,j) < limits(1,j) || new_pop(i,j) > limits(2,j)
                        new_pop(i,j) = min(new_pop(i,j),limits(2,j));
                        new_pop(i,j) = max(new_pop(i,j),limits(1,j));
                    end     
            end
        end
        % Store best value per iteration and overall
        for i=1:1:npop
            new_pop(i,4) = ff(new_pop(i,:));
        end
        
        YYY = sortrows(new_pop,-4);
        optimal_per_run(iter,:) = YYY(1,:);
            if iter == 1
                optimal_over(iter,:) = optimal_per_run(iter,:);
            elseif optimal_per_run(iter,4) >= optimal_over(iter-1,4)
                optimal_over(iter,:) = optimal_per_run(iter,:);
            else
                optimal_over(iter,:) = optimal_over(iter-1,:);
            end
        old_pop = new_pop(:,[1:3]);
    end
    final(:,run)=optimal_over(:,4);
    tabulated_results(run,:)=optimal_over(gen,:);
end

%% Display final resutls
tabulated_results(11,:)=min(tabulated_results);
tabulated_results(12,:)=max(tabulated_results);
tabulated_results(13,:)=mean(tabulated_results);
fprintf('\n            K         T1        T2       OV\n');
fprintf('Run 1 '), disp(tabulated_results(1,:)); 
fprintf('Run 2 '), disp(tabulated_results(2,:));
fprintf('Run 3 '), disp(tabulated_results(3,:)); 
fprintf('Run 4 '), disp(tabulated_results(4,:));
fprintf('Run 5 '), disp(tabulated_results(5,:)); 
fprintf('Run 6 '), disp(tabulated_results(6,:));
fprintf('Run 7 '), disp(tabulated_results(7,:)); 
fprintf('Run 8 '), disp(tabulated_results(8,:));
fprintf('Run 9 '), disp(tabulated_results(9,:)); 
fprintf('Run 10'), disp(tabulated_results(10,:));
fprintf(' Min. '), disp(tabulated_results(11,:)); 
fprintf(' Max. '), disp(tabulated_results(12,:));
fprintf(' Mean '), disp(tabulated_results(13,:));

%plot final results
iter = 300;
plot(final,'LineWidth',2), title('Real-Coded Genetic Algorithem'), ...
    xlabel('Generation','FontSize',12), ylabel('Objective function'), ...
    legend('Run 1','Run 2','Run 3','Run 4','Run 5','Run 6','Run 7','Run 8', ...
    'Run 9','Run 10'), ylim([min(min(final))-0.01 0.49])...
    , grid on;