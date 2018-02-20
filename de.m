function best = de( its )
%EATEMPLATE Summary of this function goes here
%   Detailed explanation goes here

%declare:
dimensions = 3;
population = 20;
iterations = its;
loudness = 0.95;
crossrate = 0.5;
scalefactor = 0.7;
fitness = zeros(population,1);
upper = ones(dimensions, 1)*20;
lower = ones(dimensions, 1)*-20;
fitness = zeros(population, 1);
solutions = zeros(population, dimensions);
contenders = zeros(population, dimensions);
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
best = bestfitness();
for k = 1:iterations
    for i = 1:population
        for j = 1:dimensions
            contenders(i,j) = rand() * (upper(j) - lower(j)) + lower(j);
        end
        %generate random sample of 3 members of population, not including i
        rands = randi([1 population],1,3);
        while ismember(i, rands)
            rands = randi([1 population],1,3);
        end
        %generate a random value jrand, which is a value from 1 to j.
        jrand = randi([1 dimensions],1,1);
        for j = 1:dimensions
            if rand() > crossrate || j == jrand
                %Do DE Crossover on contender
                newval = solutions(rands(1,1),j) + scalefactor * (solutions(rands(1,2),j) - solutions(rands(1,3),j));
                contenders(i,j) = bound(upper(j), lower(j), newval);
            end
        end
        tempfitness = runtrial(contenders(i,:));
        %if temp fitness is better than previous fitness, and rand is less than loudness:
        if (tempfitness < fitness(i))
            solutions(i,:) = contenders(i,:);
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
disp(solutions)
best = solutions(best,:);
end

