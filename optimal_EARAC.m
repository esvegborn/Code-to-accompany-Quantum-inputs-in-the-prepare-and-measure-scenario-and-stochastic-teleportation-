
%-------------------------------------------------------------------------%
%This function computes the average success probability of an adaptive random access code 
%assisted by a four-dimensional maximally entangled state:

% P = =\frac{1}{nY nX}\sum_{x = 1}^{nX}\sum_{y = 1}^{nY}\sum_{c=1}^{nC} \tr(A_{c \mid x} \otimes B_{x_y| y, c}\rho)

%Inputs:
% - rho: maximally entangled state
% - d = 4: local dimension
% - nX: number of Alice's settings
% - nC: number of Alice's measurement outcomes
% - nY: number of Bob's input
% - nB: number of Bob's mesurement outputs
% - A(d,d,nC,nX): Alice' local measurements
% - B(d,d,nB,nY,nC): Bob's local measurements


%Output:
% - P: the average success probability of the random access code
%-------------------------------------------------------------------------%

d=4;
nX=4^2;
nC=4;
nY=2;
nB=4;


%Creates rho (maximally entangled)
rho = zeros(1,d*d);
count = 1;
for i = 1 : d
    for j = 1 : d
        if i == j
            rho(count) = 1;
        end
        count = count + 1;
    end
end
rho = rho'*rho/d; % Default value for rho




data = load('measurements.mat');

for x = 1:nX
    for c = 1:nC
        A(:,:,c,x) = data.data{1}(:,:,c,x);
    end
end

for c = 1:nC
    for y = 1:nY
        for b = 1:nB
            B(:,:,b,y,c) = data.data{2}(:,:,b,y,c);
        end
    end
end

%Check that the measurements are valid POVMs
%-------------------------------------------------------------------------%
for x=1:nX
    s=0;
    for c=1:nC
        s=s+A(:,:,c,x);
    end
    s==eye(d);
end

for y=1:nY
    for c=1:nC
        s=0;
        for b=1:nB
            s=s+B(:,:,b,y,c);
        end
        s==eye(d);
    end
end


%Build objective function
%-------------------------------------------------------------------------%
obj=0;
for x=1:nX
    X=toSeveralBases(x-1,d*ones(1,nY))+1;
    for y=1:nY
        for c=1:nC
            obj=obj+trace(rho*Tensor(A(:,:,c,x),B(:,:,X(y),y,c)))/(nX*nY);
        end
    end
end
    

%-------------------------------------------------------------------------%
%Success probability of RAC
P = obj;


