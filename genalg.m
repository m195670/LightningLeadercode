function [MinProbBranch, MinProbDeath, MinEnergy, population] = genalg(pop, gen, n, H)
%This function uses a genetic algorithm to find the probability of
%Branching and the probability of death that will minimize the energy of
%our lightning leader

%It takes the inputs pop, which is the population size, gen which is
%the number of generations to run through, n which is the number of
%simulations to run to calculate for energy, H which is the height of the system.

[avgeng,~] = BranchingLeader(H,.5,.5,n);
test = avgeng;

%Delete the above

population = zeros(pop,3);

%This initializes the population. The first row is the Death probability,
%the second row is the Branch Probability, and the third row is the average
%energy.

for i = 1:pop
population(i, 1) = floor(rand()*1000)/1000;

population(i, 2) = floor(population(i,1)*rand()*1000)/1000;

[avgeng,~] = BranchingLeader(H, population(i,1), population(i,2), n);

population(i, 3) = abs(avgeng - test);
end

%The above for loop creates the initial population.

for i = 1:gen
    %this will continue to evaluate the population, reproduce, recombine,
    %and mutate till we reach the specified number of generations.
    
    population = sortrows(population, 3);
    
    %This evaluates and sorts the population based on the energy usage with the smallest energy being the best.
    
    reproduce = population(pop,:);
    %This is the first entry of the reproduce loop.
    for j = 1:(pop - 1)
        
        ratio = population(end,3)/population(j,3);
        
        for k = 1:ceil(ratio)
            
            reproduce(end + 1,:) = population(j,:);
        end       
    end
    %The above creates the reproduction population. These will be randomly
    %selected, recombined, and mutated in order to create a potential
    %population.
    
    potentialpop = zeros(pop, 3);
    
    for j = 1:pop 
        
        recombine = reproduce(floor(randi([1, size(reproduce,1)])),:);
        
        recombine(2,:) = reproduce(floor(randi([1, size(reproduce,1)])),:);
        
        while isequal(recombine(1,:),recombine(2,:)) == 1
            
            recombine(2,:) = reproduce(floor(randi([1,size(reproduce,1)])),:);
        end
        %The above randomly selects two members of the reproduce population
    
        crossover = de2bi(recombine(1,1)*1000,10);
            
        crossover(2,:) = de2bi(recombine(2,1)*1000,10);
       
        crossover(3,:) = zeros(1,10);
        
        crossover(3,1:4) = crossover(1,1:4);
        
        crossover(3, 5:10) = crossover(2,5:10);
        
        crossover(4,:) = de2bi(recombine(1,2)*1000,10);
            
        crossover(5,:) = de2bi(recombine(2,2)*1000,10);
       
        crossover(6,:) = zeros(1,10);
        
        crossover(6,1:4) = crossover(4,1:4);
        
        crossover(6, 5:10) = crossover(5,5:10);
        
        %The above crossesover the first 4 bits for the death prob and
        %branch prob.
        
        for k = 1:10
           
            x = rand();
            
            if x <= .03
                
               crossover(3,k) = mod(crossover(3,k) + 1, 2); 
            end
            
        end
        
        for k = 1:10
           
            x = rand();
            
            if x <= .03
                
               crossover(6,k) = mod(crossover(6,k) + 1, 2); 
            end
            
        end
        %The above mutates each bit with probability 3%
        
        potentialpop(j,1) = bi2de(crossover(3,:))/1000;
        
        potentialpop(j,2) = bi2de(crossover(6,:))/1000;
        
        if potentialpop(j,1) + potentialpop(j,2) > 1
            
            potentialpop(j,3) = inf;
        
        elseif sum(ismember(population(:,1:2), potentialpop(j, 1:2))) > 0
            potentialpop(j,3) = inf;
        else 
        [avgeng,~] = BranchingLeader(H, potentialpop(j,1), potentialpop(j,2), n);
        
        potentialpop(j,3) = abs(avgeng - test);
        
        %The above calculates the fitness of the potential population
        end
    end

    for j = 1:size(potentialpop,1)
       
        if potentialpop(j,3) < population(end, 3)
            
            population(end, :) = potentialpop(j,:);
            
        end
        
        population = sortrows(population, 3);
        
    end
    
    
end
MinProbBranch = population(1,2);

MinProbDeath = population(1,1);

MinEnergy = population(1,3);


end
