function [MinProbDeath, MinProbBranch, MinEnergy, path] = OptimumEng(H, Death, Branch)
format longG
if Death + Branch > 1
    
    fprintf("inputs are outside of the bounds\n")
    
    MinProbDeath = "undefined";
    
    MinProbBranch = "undefined";
    
    MinEnergy = "undefined";
    
    path = "undefined";
    
    return
end

[testing, ~] = BranchingLeader(H, .8, .02, 4);
%delete^
path = [Death, Branch];
%creates the path that the code is taking to optimize
n = 4;
%change to 1000000^
%sets the initial simulations to run
[avgeng, Energy, ~] = BranchingLeader(H, Death, Branch, n);
%runs the code to get the average energy
B = 0;

    for j = 1:n
        B = B + (Energy(j) - avgeng)^2;
    end
    
S = sqrt((1/(n-1))*B);

A = (1.96*S)/(sqrt(n))

%The above creates the confidence interval

while A > 100000000
    %change to 1^
   n = 2*n;
   
   [avgeng, Energy, ~] = BranchingLeader(H, Death, Branch, n);
   
B = 0; 

   for j = 1:n
        B = B + (Energy(j) - avgeng)^2;
   end
    
S = sqrt((1/(n-1))*B);

A = (1.96*S)/(sqrt(n));
end
%The above verifies that the confidence interaval is less than 1, if it is not it will redo the simulation with double the runs.    
NewEnergy = abs(avgeng - testing);
%change to avgeng^
%changes NewEnergy to the average energy
patheng = NewEnergy;
%updates the patheng to record the energy
OldEnergy = NewEnergy;
%sets the old energy to the new energy
input = zeros(4,2);
%this creates the input matrix, which should indicate all of the new inputs
input(1:2, 2) = Branch;
%simply sets two of the rows to have the same Branch probability
input(3:4, 1) = Death;
%sets the other two rows to have the same death probability
input(1,1) = Death + .5;
%sets the first input to check [Branch, Death + .5]
input(2,1) = Death - .5;
%sets the second input to check [Branch, Death - .5]
input(3,2) = Branch + .5;
%sets the Third input to check [Branch + .5, Death]
input(4,2) = Branch - .5;
%sets the fourth input to check [Branch + .5, Death]
TestedEnergy = zeros(1,4);
%creates a row vector to store the tested energy for all of the inputs
for i = 1:4
   if input(i,1) + input(i,2) > 1
       
       TestedEnergy(i) = inf;
       
       %this says that if the probability of branching and the probability
       %of death are greater than 1
   elseif (input(i,1) >= 0) && (input(i,1) <= 1) && (input(i,2) >= 0) && (input(i,2) <= 1)
       %states that we are only testing the new inputs if they are within
       %the probability space
       n = 4;
%change to 1000000^
%sets the initial number of simulations to run
       [avgeng, Energy, ~] = BranchingLeader(H, input(i,1) , input(i,2) , n);
%Runs the simulation with the desired inputs
B = 0;

       for j = 1:n
           B = B + (Energy(j) - avgeng)^2;
       end
    
       S = sqrt((1/(n-1))*B);

       A = (1.96*S)/(sqrt(n));
%sets the initial confidence interaval
       while A > 100000000
    %change to 1^
           n = 2*n;
   
           [avgeng, Energy, ~] = BranchingLeader(H, input(i,1), input(i,2), n);

           B = 0;
           
           for j = 1:n
               B = B + (Energy(j) - avgeng)^2;
           end
    
       S = sqrt((1/(n-1))*B);

       A = (1.96*S)/(sqrt(n));
       
       end
%makes sure that the confidence interval is less than 1 otherwise reruns the test with double the runs      
       TestedEnergy(i) = abs(avgeng - testing);
       %change to avgeng^
       %stores the energy in our tested energy array
   else
       
       TestedEnergy(i) = inf;
   %if the selected probabilities are not in the space we are setting the
   %energy to infinity so that the path isn't chosen
   end
   
   
end
[Minimum, index] = min(TestedEnergy);
%finding the minimum and where the minimum index is in our energy vector
NewEnergy = Minimum;
%setting new energy equal to the minimum tested energy
if NewEnergy < OldEnergy
patheng(2,1) = NewEnergy;
%updates the energy path with the new energy
path(2,:) = input(index,:);
%updates the path with the next step
z = 3;
else 
z = 2;
end

%Simply a counter to update the path indexes
while NewEnergy < OldEnergy
    %stating that we are only running through this loop while the new
    %energy is less than the old energy
OldEnergy = NewEnergy;
   %this updates the old energy to equal the new energy so that in the next
   %iteration we can see if the new energy is getting smaller or not.
input = zeros(4,2);
%rewrites the input matrix
input(1:2, 2) = path(z-1,2);
%updates so that the Branch probabilities are the same as the last step
input(3:4, 1) = path(z-1,1);
%updates so that the Death probabilities are the same as the last step
input(1,1) = path(z-1, 1) + .5;

input(2,1) = path(z-1, 1) - .5;

input(3,2) = path(z-1, 2) + .5;

input(4,2) = path(z-1, 2) - .5;
%updates the new inputs with the appropriate steps to take

TestedEnergy = zeros(1,4);
%rewrites the energy that we will calculate

for i = 1:4
   if 1 == any(ismember(path,input(i,:), 'rows'))
       
       TestedEnergy(i) = inf;
       %if the input is already a member of the path, we will not even test
       %it, this should save a lot of time since if it is a member of the
       %path it shouldn't be lower than the current energy
       
   elseif input(i,1) + input(i,2) > 1
       
       TestedEnergy(i) = inf;
       
       %this says that if the probability of branching and the probability
       %of death are greater than 1
   elseif (input(i,1) > 0) && (input(i,1) < 1) && (input(i,2) > 0) && (input(i,2) < 1)
       %verifies that the input is inside of the probability space
       n = 4;
%change to 1000000^
       [avgeng, Energy, ~] = BranchingLeader(H, input(i,1) , input(i,2) , n);

       B = 0;
       
       for j = 1:n
           B = B + (Energy(j) - avgeng)^2;
       end
    
       S = sqrt((1/(n-1))*B);

       A = (1.96*S)/(sqrt(n));
%sets the confidence interval
       while A > 100000000
    %change to 1^
           n = 2*n;
   
           [avgeng, Energy, ~] = BranchingLeader(H, input(i,1), input(i,2), n);
   
           B = 0;
           
           for j = 1:n
               B = B + (Energy(j) - avgeng)^2;
           end
    
       S = sqrt((1/(n-1))*B);

       A = (1.96*S)/(sqrt(n));
       %makes sure that the confidence level is below 1
       end
       
       TestedEnergy(i) = abs(avgeng - testing);
       %updates the tested energy
   else
       
       TestedEnergy(i) = inf;
   %sets the tested energy to infinity if the path is outside of the
   %bounds.
   end
   
   
end
[Minimum, index] = min(TestedEnergy);
%finds the minimum energy
NewEnergy = Minimum;
%sets the new energy to the minimum
if NewEnergy < OldEnergy
path(z,:) = input(index,:);
%updates the path
patheng(z,1) = NewEnergy;
%updates the energypath
z = z + 1;
%updates the index counter
end

end
searchradius = [.25, .1, .05, .01];
%This is the array of the different search radius.
for k = 1:4
input = zeros(4,2);
%rewrites the input matrix
input(1:2, 2) = path(z-1,2);
%updates so that the Branch probabilities are the same as the last step
input(3:4, 1) = path(z-1,1);
%updates so that the Death probabilities are the same as the last step
input(1,1) = path(z-1, 1) + searchradius(k);

input(2,1) = path(z-1, 1) - searchradius(k);

input(3,2) = path(z-1, 2) + searchradius(k);

input(4,2) = path(z-1, 2) - searchradius(k);
%updates the new inputs with the appropriate steps to take
TestedEnergy = zeros(1,4);
%rewrites the energy that we will calculate
for i = 1:4
   if 1 == any(ismember(path,input(i,:), 'rows'))
       
       TestedEnergy(i) = inf;
       %if the input is already a member of the path, we will not even test
       %it, this should save a lot of time since if it is a member of the
       %path it shouldn't be lower than the current energy
       
   elseif input(i,1) + input(i,2) > 1
       
       TestedEnergy(i) = inf;
       
       %this says that if the probability of branching and the probability
       %of death are greater than 1   
       
   elseif (input(i,1) >= 0) && (input(i,1) <= 1) && (input(i,2) >= 0) && (input(i,2) <= 1)
       %verifies that the input is inside of the probability space
       n = 4;
%change to 1000000^
       [avgeng, Energy, ~] = BranchingLeader(H, input(i,1) , input(i,2) , n);

       B = 0;
       
       for j = 1:n
           B = B + (Energy(j) - avgeng)^2;
       end
    
       S = sqrt((1/(n-1))*B);

       A = (1.96*S)/(sqrt(n));
%sets the confidence interval
       while A > 100000000
    %change to 1^
           n = 2*n;
   
           [avgeng, Energy, ~] = BranchingLeader(H, input(i,1), input(i,2), n);
   
           B = 0;
           
           for j = 1:n
               B = B + (Energy(j) - avgeng)^2;
           end
    
       S = sqrt((1/(n-1))*B);

       A = (1.96*S)/(sqrt(n));
       %makes sure that the confidence level is below 1
       end
       
       TestedEnergy(i) = abs(avgeng - testing);
       %updates the tested energy
   else
       
       TestedEnergy(i) = inf;
   %sets the tested energy to infinity if the path is outside of the
   %bounds.
   end
   
   
end
[Minimum, index] = min(TestedEnergy);
%finds the minimum energy
NewEnergy = Minimum;
%sets the new energy to the minimum
if NewEnergy < OldEnergy
path(z,:) = input(index,:);
%updates the path
patheng(z,1) = NewEnergy;
%updates the energypath
z = z + 1;
%updates the index counter
end
while NewEnergy < OldEnergy
    %stating that we are only running through this loop while the new
    %energy is less than the old energy
OldEnergy = NewEnergy;
   %this updates the old energy to equal the new energy so that in the next
   %iteration we can see if the new energy is getting smaller or not.
input = zeros(4,2);
%rewrites the input matrix
input(1:2, 2) = path(z-1,2);
%updates so that the Branch probabilities are the same as the last step
input(3:4, 1) = path(z-1,1);
%updates so that the Death probabilities are the same as the last step
input(1,1) = path(z-1, 1) + searchradius(k);

input(2,1) = path(z-1, 1) - searchradius(k);

input(3,2) = path(z-1, 2) + searchradius(k);

input(4,2) = path(z-1, 2) - searchradius(k);
%updates the new inputs with the appropriate steps to take
TestedEnergy = zeros(1,4);
%rewrites the energy that we will calculate
for i = 1:4
   if 1 == any(ismember(path,input(i,:), 'rows'))
       
       TestedEnergy(i) = inf;
       %if the input is already a member of the path, we will not even test
       %it, this should save a lot of time since if it is a member of the
       %path it shouldn't be lower than the current energy
   elseif input(i,1) + input(i,2) > 1
       
       TestedEnergy(i) = inf;
       
       %this says that if the probability of branching and the probability
       %of death are greater than 1
   elseif (input(i,1) > 0) && (input(i,1) < 1) && (input(i,2) > 0) && (input(i,2) < 1)
       %verifies that the input is inside of the probability space
       n = 4;
%change to 1000000^
       [avgeng, Energy, ~] = BranchingLeader(H, input(i,1) , input(i,2) , n);

       B = 0;
       
       for j = 1:n
           B = B + (Energy(j) - avgeng)^2;
       end
    
       S = sqrt((1/(n-1))*B);

       A = (1.96*S)/(sqrt(n));
%sets the confidence interval
       while A > 100000000
    %change to 1^
           n = 2*n;
   
           [avgeng, Energy, ~] = BranchingLeader(H, input(i,1), input(i,2), n);
   
           B = 0;
           
           for j = 1:n
               B = B + (Energy(j) - avgeng)^2;
           end
    
       S = sqrt((1/(n-1))*B);

       A = (1.96*S)/(sqrt(n));
       %makes sure that the confidence level is below 1
       end
       
       TestedEnergy(i) = abs(avgeng - testing);
       %updates the tested energy
   else
       
       TestedEnergy(i) = inf;
   %sets the tested energy to infinity if the path is outside of the
   %bounds.
   end
   
   
end
[Minimum, index] = min(TestedEnergy);
%finds the minimum energy
NewEnergy = Minimum;
%sets the new energy to the minimum
if NewEnergy < OldEnergy
path(z,:) = input(index,:);
%updates the path
patheng(z,1) = NewEnergy;
%updates the energypath
z = z + 1;
%updates the index counter
end
end
end


MinProbDeath = path(end, 1);
%shows the optimum Probability of Death
MinProbBranch = path(end, 2);
%shows the optimum probability of Branching
MinEnergy = OldEnergy;
%shows what the minimum energy calculated was

path = [path, patheng];
%changes path to output the energy as well as teh different probabilities
%tested.
