function best = ba( its )
%EATEMPLATE Summary of this function goes here
%   Detailed explanation goes here

%declare:
dimensions = 3;
population = 20;
iterations = its;
loudness = 0.95;
pulserate = 0.1;
qmin = 0;
qmax = 2;
fitness = zeros(population,1);
upper = ones(dimensions, 1)*20;
lower = ones(dimensions, 1)*-20;
fitness = zeros(population, 1);
solutions = zeros(population, dimensions);
velocity = zeros(population, dimensions);
tempsol = zeros(population, dimensions);
for i = 1:population
    for j = 1:dimensions
        solutions(i,j) = rand() * (upper(j) - lower(j)) + lower(j);
    end
end
dir = char(datetime('now','Format','yyyy-MM-dd''T''HHmmss'));
mkdir(dir);

function fitness = runtrial(bat)
    fitness = sqrt((bat(1))^2) + sqrt((bat(2))^2) + sqrt((bat(3))^2);
end

function value = bound(upper, lower, val)
    if val > upper
        value = upper;
    elseif val < lower
        value = lower;        
    else
        value = val;
    end
end

function updatefitnesses()
    for m = 1:population
        fitness(m) = runtrial(solutions(m,:));
    end
end

function best = bestfitness()
    b = 1;
    for n = 1:population
        if fitness(n) < fitness(b)
            b = n;
        end
    end
    best = b;
end

updatefitnesses()
for k = 1:iterations
    q = rand(population,1) * (qmax - qmin) + qmin;
    best = bestfitness();
    %solutions(1,:)
    %disp('Start of loop')
    for i = 1:population
        velocity(i,:) = velocity(i,:) + (solutions(i,:) - solutions(best,:))*q(i);
        disp(velocity)
        tempsol(i,:) = solutions(i,:) + velocity(i,:);
            %disp(tempsol(i,j))
        for j = 1:dimensions
            tempsol(i,j) = bound(upper(j), lower(j), tempsol(i,j));
        end
            %if i == 1
            %    solutions(1,:)
            %    disp('post tempsol update')    
            %end
        % if rand greater than pulse rate:
        if rand() > pulserate
            for j = 1:dimensions
                newval = stblrnd(1,1,4,8)*0.001 * tempsol(i,j);
                tempsol(i,j) = bound(upper(j), lower(j), newval);
                %if i == 1
                %    solutions(1,:)
                %    disp('The pulse happened')
                %end
            end
        end
        tempfitness = runtrial(tempsol(i,:));
        %if temp fitness is better than previous fitness, and rand is less than loudness:
        if (tempfitness < fitness(i)) && (rand() < loudness)
            %update the actual position with the temp position
            solutions(i,:) = tempsol(i,:);
            %if i == 1
            %    solutions(1,:)
            %    disp('Solutions was updated')
            %end
            fitness(i) = tempfitness;
        end
        %if it's better than best fitness, update best fitness.
        if tempfitness < fitness(best)
            best = i;
        end
        %solutions(1,:)
        %disp('End of loop')
    end
    figure('visible','off')
    x = solutions(:,1);
    y = solutions(:,2);
    z = solutions(:,3);
    scaledfit = fitness - min(fitness);
    scaledfit = scaledfit / max(scaledfit);
    graph = scatter3(x,y,z,[],scaledfit,'filled');
    colorbar
    file = int2str(k);
    name = '.png';
    filename = strcat(file,name);
    fname = strcat('C:\Users\Ben\Documents\MATLAB\',dir);
    saveas(graph, fullfile(fname, filename));
end
best = solutions(best,:);
end

