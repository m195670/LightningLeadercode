function [B, avgenergy, Energy, avgsteps, Steps, Exp_p, Var_p] = BranchingLeader(H,D,Br,n)
tic
%H = starting Height of energy

%D = Probability of Death (to nearest .001)

%Br = Probability of Branching (to nearest .001)

%n = number of times to run simulation

pd = D*1000;

%This is the probability of death range, used in the random variate
%generator later on

pbr = pd + Br*1000;

%This is the probability of branching, used in the random variate generator
%later on

Steps = zeros(1,n);
   
    %Steps is the total number of steps that was taken before the leader
    %hit the ground or reached our maximum output.
    
   
Energy = zeros(1,n);
    
% Energy is the "Energy" expended with the leader as a whole
    


for i = 1:n
    
    %This will run the code n number of times and store the values
    
    E = 0;

    %E is the energy used
    
    T = 0;

    %T is the step number

    B = zeros(1,2);
    
    %B is the vector that the leaders will take. 
    
    B(1,2) = H;
    
    %This states initializes B for the Height. 
    
    Leaders = zeros(1,1000);
    
    %Leaders is a binary matrix that indicates when a leader is on.
    
    Leaders(1,1) = 1;
    
    %This indicates that the first leader is 'alive'
    
    Z = randi([1, 3]);
    
    %Random variate generator, should output {1,2,3} with equal probability
  
    if Z == 1
    
        B(2,:) = B(1,:) + [0, -1];
   
    elseif Z == 2
    
        B(2,:) = B(1,:) + [1, 0];
    
    else
        
        B(2,:) = B(1,:) + [-1, 0];
    
    end
    
    %this takes the output of a random variate generator and decides which
    %direction to go
    
    T = T + 1;
    
    %This accounts for the step that was just taken
    
    A = B(2,2);
    
    %This will indicate the lowest Height. 
    
    w = 3;
    
    %This is used as an index
    
    E = E + 1;
    
    %Energy expended in the last step
    
    Leaders(1,2) = 1;
    
        while (A ~= 0) 
            

            if T >= 999
                
                break
                
            end
            
            %This will break the while loop if we go over 1000 steps.
            
            
            for j = 1:size(Leaders,1)
                
            %This will cycle through all of the current leaders that have
            %been created
             
            B_D = randi([0,1000]);
            
                %This is a random variate generator for determining Branching
                %or Death. 
                
                if Leaders(j,w-1) == 0 
                
                    Leaders(j,w) = 0;
            
                    %This specifies that if the leader has died, it will not come
                    %back.
                
                    B(w,2*j-1:2*j) = inf;
                    
                    %This sets the position infinity.
            
                elseif sum(Leaders(1:j-1,w)) + sum(Leaders(j:end,w-1)) == 1
                    
                    BR = Br*1000;
                    S = 1000 - pbr;
                    Max = S + BR;
                    B_S = randi([0, floor(abs(Max))]);
                    if (0 <= B_S) && ( B_S <= BR)
                        Last_Step = B(w-1,2*j-1:2*j)-B(w-2,2*j-1:2*j);
            
                    %This says that the last step is what is stored in last entry
                    %of our B matrix subracted from our 2nd to last entry. Should
                    %be an entry [1,0], [-1,0], or [0,-1].
                    
                    Leaders(j,w) = 1;
                    
                    %This tells the matrix that the leader is still 'alive'
                    
                    Leaders(end+1,1:w) = 0;
                    
                    %This creates another row for leaders to indicate the
                    %branch
                    
                    Leaders(end, w) = 1;
                    
                    %This tells Leaders that the created branch is alive.
            
                    B(1:w, end+1:end+2) = 0;
                    
                    %This initializes the first steps to zero
                    
                    B(w - 1, end-1:end) = B(w-1, 2*j-1:2*j);
                    
                    Z = randi([1,3]);
                    
                    %Random integer between 1 and 3, determines which step
                    %is taken
                    
                    if Z==1
                        
                        B(w, 2*j-1:2*j) = B(w-1,2*j-1:2*j) + Last_Step;
                        
                        B(w, end-1:end) = B(w-1, 2*j-1:2*j) + [Last_Step(2), -Last_Step(1)];
                        
                    elseif Z==2
                        
                        B(w, 2*j-1:2*j) = B(w-1,2*j-1:2*j) + Last_Step;
                        
                        B(w, end-1:end) = B(w-1, 2*j-1:2*j) - [Last_Step(2),-Last_Step(1)];
                    else
                        B(w, 2*j-1:2*j) = B(w-1, 2*j-1:2*j) - [Last_Step(2), -Last_Step(1)];
                        
                        B(w, end-1:end) = B(w-1, 2*j-1:2*j) + [Last_Step(2), -Last_Step(1)];
                    end
                    
                    %This goes through the options for the next step to
                    %take
                    
                    E = E + 2;
                    
                    %Adds the amount of energy for two leaders steping
                    %once.
                    else
                        Last_Step = B(w-1,2*j-1:2*j)-B(w-2,2*j-1:2*j);
                    Z = randi([1, 3]);
                       if Z==1
                           
                           B(w, 2*j-1:2*j) = B(w-1, 2*j-1:2*j) + Last_Step;
                           
                        elseif Z==2
                            
                            B(w, 2*j-1:2*j) = B(w-1, 2*j-1:2*j) + [Last_Step(2),-Last_Step(1)];
    
                       else
                           
                            B(w, 2*j-1:2*j) = B(w-1,2*j-1:2*j) - [Last_Step(2),-Last_Step(1)];
                        
                       end 
                       
                       E = E + 1;
                       
                       %Energy used by one leader stepping once
                       
                       Leaders(j,w) = 1;
                    end   
                    
                elseif (0 <= B_D) && ( B_D <= pd)
                
                    Leaders(j,w) = 0;
                    
                    B(w,2*j-1:2*j) = inf;
                
                elseif (B_D > pd) && (B_D <= pbr)
                
                    Last_Step = B(w-1,2*j-1:2*j)-B(w-2,2*j-1:2*j);
            
                    %This says that the last step is what is stored in last entry
                    %of our B matrix subracted from our 2nd to last entry. Should
                    %be an entry [1,0], [-1,0], or [0,-1].
                    
                    Leaders(j,w) = 1;
                    
                    %This tells the matrix that the leader is still 'alive'
                    
                    Leaders(end+1,1:w) = 0;
                    
                    %This creates another row for leaders to indicate the
                    %branch
                    
                    Leaders(end, w) = 1;
                    
                    %This tells Leaders that the created branch is alive.
            
                    B(1:w, end+1:end+2) = 0;
                    
                    %This initializes the first steps to zero
                    
                    B(w - 1, end-1:end) = B(w-1, 2*j-1:2*j);
                    
                    Z = randi([1,3]);
                    
                    %Random integer between 1 and 3, determines which step
                    %is taken
                    
                    if Z==1
                        
                        B(w, 2*j-1:2*j) = B(w-1,2*j-1:2*j) + Last_Step;
                        
                        B(w, end-1:end) = B(w-1, 2*j-1:2*j) + [Last_Step(2), -Last_Step(1)];
                        
                    elseif Z==2
                        
                        B(w, 2*j-1:2*j) = B(w-1,2*j-1:2*j) + Last_Step;
                        
                        B(w, end-1:end) = B(w-1, 2*j-1:2*j) - [Last_Step(2),-Last_Step(1)];
                    else
                        B(w, 2*j-1:2*j) = B(w-1, 2*j-1:2*j) - [Last_Step(2), -Last_Step(1)];
                        
                        B(w, end-1:end) = B(w-1, 2*j-1:2*j) + [Last_Step(2), -Last_Step(1)];
                    end
                    
                    %This goes through the options for the next step to
                    %take
                    
                    E = E + 2;
                    
                    %Adds the amount of energy for two leaders steping
                    %once.

                else
                    Last_Step = B(w-1,2*j-1:2*j)-B(w-2,2*j-1:2*j);
                    Z = randi([1, 3]);
                       if Z==1
                           
                           B(w, 2*j-1:2*j) = B(w-1, 2*j-1:2*j) + Last_Step;
                           
                        elseif Z==2
                            
                            B(w, 2*j-1:2*j) = B(w-1, 2*j-1:2*j) + [Last_Step(2),-Last_Step(1)];
    
                       else
                           
                            B(w, 2*j-1:2*j) = B(w-1,2*j-1:2*j) - [Last_Step(2),-Last_Step(1)];
                        
                       end 
                       
                       E = E + 1;
                       
                       %Energy used by one leader stepping once
                       
                       Leaders(j,w) = 1;

                  
                end
            
            end
            
            A = min(B(w, 2:2:end));
            
            w = w + 1;
            
            %This creates the indicies. Its also the amount of steps taken
            %so far
                                  
            T = T + 1;
        end
        Energy(i) = E;
        
        %This sets the Energy for each level of the for loop
        
        Steps(i) = T; 
        
        %This sets the steps for each level of the for loop 
        
        if A == 0
            
        [~,I] = min(B(end,2:2:end));
        
        V(i) = B(end, 2*I - 1);
        
        else
          
        V(i) = inf;
        
        
        end
        
        %Stores where the Leader hits only if the leader hits the ground
end

avgenergy = sum(Energy)/n;

%average energy for the entire thing.

avgsteps = sum(Steps)/n;

%average steps for the entire thing.

V = V(V~=inf);

%This take out all of the values in which the leader did not hit the
%ground.

%Probground = size(V,2)/n;

p_bins = min(V):max(V);

%makes bins that contain all integers between the minimum and maximum

%[count_V] = hist(V,p_bins);

%This counts all of the points in the histogram that correspond to each of
%the integers in p_bins

%pdf_x = count_V/n;

%This is the pmf

%plot(p_bins ,pdf_x);

%Just plots the probabilities

%Exp_p = sum(p_bins .* pdf_x);

%Calculates the expected value

%Var_p = sum((p_bins.^2).* pdf_x) - (Exp_p)^2;

%Calculates the expected value
toc
end
