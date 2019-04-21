function [Probability] = nLightningLeaders(H,D)
%Returns the probability that a spot is hit by at least one lightning
%leader

%H is the height, D is a row vector of distances between starting
%points

N = H-1;

%N is the number of steps

TD = sum(D);

%Total distance between the starting positions for the Lightning Leader

n = length(D)+1;

%length of input vector

P = zeros(n, H + TD);

%this creates a matrix that has each alpha leader as a row w/ row size
%being H plus the total distance added by the space between start points

for i = 1:(N+1)
    P(1, i) = nchoosek(N,i-1)*(1/2)^(N);
end

%This creates the first row of probabilities

for Z = 2:n
   for i = 1:(N+1)
        P(Z, i + sum(D(1:Z-1))) = nchoosek(N,i-1)*(1/2)^(N);
   end
end

%This creates the rest of the rows of probabilities

Probability = sum(P(1:n,1:N+1+TD));

w = 2;

while w <= n
    
    A = combnk(1:n,w);
    
    D = size(A);
    
    C = zeros(1, N + 1 + TD);
    
    for i = 1:D(1)
        
        B = ones(1, N + 1 + TD);
        
        for j = 1:D(2)
            
            B = B.*P(A(i,j),:);
            
        end
        
        C = C + B;
        
    end 
    
Probability = Probability + (-1)^(w-1).*C;

w = w + 1;

end

%The above applies the law of total probability

end