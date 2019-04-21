function [Allavgeng, AllEnergy, Allavgsteps, AllSteps] = BruteForce(Branch, m)
AllEnergy = zeros(21 , m + 1);
AllEnergy(:,1) = 0:.01:.2;
Allavgeng = zeros(21,3);
Allavgeng(:,1) = 0:.01:.2;
Allavgsteps = zeros(21,3);
Allavgsteps(:,1) = 0:.01:.2;
AllSteps = zeros(21, m + 1);
AllSteps(:,1) = 0:.01:.2;
n = 1;

for i = 0:.01:.2
    B = 0;
    C = 0;
   [Avgeng,Energy,Avgsteps, Steps] = BranchingLeader(10,i,Branch, m);
   Allavgeng(n,2) = Avgeng;
   Allavgsteps(n,2) = Avgsteps;
    for j = 1:m
        B = B + (Energy(j) - Avgeng)^2;
    end
    S = sqrt((1/(m-1))*B);
    A = (1.96*S)/(sqrt(m));
    Allavgeng(n,3) = A;
    for j = 1:m
        C = C + (Steps(j) - Avgsteps)^2;
    end
    D = sqrt((1/(m-1))*C);
    E = (1.96*D)/(sqrt(m));
    Allavgsteps(n,3) = E;
    AllEnergy(n,2:end) = Energy;
    AllSteps(n, 2:end) = Steps;
    n = n + 1;
end
Allavgeng = array2table(Allavgeng, 'VariableNames',{'DeathProbabilities','AverageEnergy', 'ConfidenceInterval'});
Allavgsteps = array2table(Allavgsteps, 'VariableNames',{'DeathProbabilities','AverageSteps', 'ConfidenceInterval'});
name = strcat("Average Energy for ", num2str(Branch), " Prob Branching");
save(name)