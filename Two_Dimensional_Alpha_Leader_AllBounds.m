function [Exp_x, Var_x] = Two_Dimensional_Alpha_Leader_AllBounds(H, n)
%Creates a row vector to capture where on the leader hits the ground
V = zeros(1,n);

%Runs through the code n number of times
for i = 1:n
%This creates B which is the path the lightning takes
B = zeros(1,2);

%Sets the leader to start at (0,h)
B(1,2) = H;

%Random variate generator, should output {1,2,3} with equal probability
Z = ceil(3*rand());

%this takes the output of the random variate generator and decides which
%direction to go
if Z == 1
    B(2,:) = B(1,:) + [0, -1];
elseif Z == 2
    B(2,:) = B(1,:) + [1, 0];
else 
    B(2,:) = B(1,:) + [-1, 0];
end

%This is just an indice that updates after each iteration of the while loop
w = 3;
A = 1;
%Should create the path to ground (Runs through until the leader hits the
%ground
while (A ~= 0)
    %Specifies which direction to go, which is dependent on the previous
    %step
    
    if w >1000
            
        B(3:end,:) = [];
        
        w = 3;
        
    end
    
    Last_Step = B(w-1,:)-B(w-2,:);
    if B(w-1,2) == H
        
        if isequal(Last_Step,[0,1])
            Z = randi([1,2]);
            if Z == 1
                B(w,:) = B(w-1,:) + [Last_Step(2),-Last_Step(1)];
            else 
                B(w,:) = B(w-1,:) - [Last_Step(2),-Last_Step(1)];
            end
        else 
            Z = randi([1,2]);
            if Z == 1
                 B(w,:) = B(w-1,:) + Last_Step;
            else
                 B(w,:) = B(w-1,:) + [0, -1];
            end
        end
    elseif abs(B(w-1,1)) == 2*H
        if isequal(Last_Step, [1,0]) || isequal(Last_Step,[-1,0])
            Z = randi([1,2]);
            if Z == 1
                B(w,:) = B(w-1,:) + [Last_Step(2),-Last_Step(1)];
            else 
                B(w,:) = B(w-1,:) - [Last_Step(2),-Last_Step(1)];
            end
        elseif B(w-1,1) == 2*H 
            Z = randi([1,2]);
            if Z == 1
                 B(w,:) = B(w-1,:) + Last_Step;
            else
                 B(w,:) = B(w-1,:) + [-1, 0];
            end
        else
                        Z = randi([1,2]);
            if Z == 1
                 B(w,:) = B(w-1,:) + Last_Step;
            else
                 B(w,:) = B(w-1,:) + [1, 0];
            end
        end
    else
    Z = randi([1,3]);
    if Z==1
        B(w,:) = B(w-1,:) + Last_Step;
    elseif Z==2
        B(w,:) = B(w-1,:) + [Last_Step(2),-Last_Step(1)];
    else
        B(w,:) = B(w-1,:) - [Last_Step(2),-Last_Step(1)];
    end
    end
  A = B(w,2);
%Just for the indices
   w = w + 1; 
end

% This records the place the leader hit the ground
V(i) = B(end, 1);
%hold on
%plot(B(:,1),B(:,2))

end

%makes bins that contain all integers between the minimum and maximum
p_bins = min(V):max(V);

%This counts all of the points in the histogram that correspond to each of
%the integers in p_bins
[count_V] = hist(V,p_bins);

%This is the pdf 
pdf_x = count_V/n;

%Just plots the probabilities
plot(p_bins ,pdf_x);

%Calculates the expected value
Exp_x = sum(p_bins .* pdf_x);

Var_x = sum((p_bins.^2).* pdf_x) - (Exp_x)^2;

end