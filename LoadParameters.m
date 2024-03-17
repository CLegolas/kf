% Parameter setting
mag_min  = 10; 
mag_max  = 30;
angle_min = 0; 
angle_max = 2*pi;

m1 = sqrt(4);
for i=1:4
    if mod(i,m1)==0
        pos(i,1)=500;
        pos(i,2)=(floor(i/m1)-1)/(m1-1)*500;
    else
        pos(i,1)=(mod(i,m1)-1)/(m1-1)*500;
        pos(i,2)=floor(i/m1)/(m1-1)*500;
    end
end
% network connectivity
degree = 2; 
E = getConnectivity(4,degree); 
% for ii=1:4
%     if ii<=4-1
%         E(ii,ii+1)=1;
%     end
%     if ii>=2
%         E(ii,ii-1)=1;
%     end
% end
%FILTER parameter
R = 100*eye(2,2); 
H = [1,0,0,0; 0,1,0,0]; 
Rinv = inv(R);
meanc =[0;0]; 
Phi = [1, 0, 1, 0;  0, 1, 0, 1;   0, 0, 1, 0;  0, 0, 0, 1 ]; 
Q = diag([10, 10, 1 , 1]); 
Tr(:,:,1) = zeros(4,4) ;
Tr(1,1,1) = 1 ;
Tr(3,3,1) = 1 ;
Tr(:,:,2) = zeros(4,4) ;
Tr(2,2,2) = 1 ;
Tr(4,4,2) = 1 ;
%%
% E=zeros(4,4);
% for ii=1:4
%     if ii<=4-1
%         E(ii,ii+1)=1;
%     end
%     if ii>=2
%         E(ii,ii-1)=1;
%     end
%     
% end

%% Initialization
xa = zeros(4,100); 
while 1
 flag_restart = 0;
 sp =500*rand;
  dir = angle_min+(angle_max-angle_min)*rand;
  x = 500 * rand;
  y = 500 * rand;
xa(:,1) = [x , y , sp*cos(dir) , sp*sin(dir)  ]';
    for t = 1:100
        if t>1
        xa(:,t) = Phi*xa(:,t-1)+mvnrnd(zeros(4,1),Q)';
         if xa(1,t)> 500 || xa(1,t) < 0 || xa(2,t) > 500 || xa(2,t) < 0
         flag_restart = 1;
         break;
  end
        end
     end 
    if flag_restart
     continue;
      elseif norm( xa(1:2,1) - xa(1:2,100) ) < 200 %
      continue;
    else
     break;
end
end

P = cell(4,1);
eta = cell(4,1);

for i=1:4
    P{i} = diag([100,100,10,10]);
    eta{i} = mvnrnd(zeros(4,1),P{i})';
end

zt = cell(100,1);
zCountt = cell(100,1);
meanc = [0;0];
for EE=1:100
    z = cell(4,1);
    zCount = zeros(4,1);
    for i = 1:4
        ztemp = H*xa(:,EE);
        z{i} = ztemp+mvnrnd(meanc,R)';
        zCount(i) = zCount(i) + 1; 
    end
    zt{EE} = z;
    zCountt{EE} = zCount;
end