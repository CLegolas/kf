% Code Author: Cao Ruihao - 465964222@qq.com
% no permission necessary for non-commercial use
% Date: 3/13/2034

%% Grid example¡ª¡ªMultiple sensors are placed according to a grid
clc;
clear;
close all

commu_iteration = 22/2;
time_rgicf  = zeros(1,commu_iteration);
time_ggicf  = zeros(1,commu_iteration);
time_sggicf = zeros(1,commu_iteration);
time_icf     = zeros(1,commu_iteration);
time_icf_P_U = zeros(1,commu_iteration);

    run LoadParameters
    for nn = (1:commu_iteration)*1
        K=2*nn
        threshold = 0.5;           
        eps_constant = .65/(max(sum(E))); 
        eps=eps_constant;

        rgicf  =  RGICF(eta,xa,P); 
        ggicf  =  GGICF(eta,xa,P); 
        sggicf = SGICF(eta,xa,P); 
        icf     =     ICF(eta,xa,P); 
        icf_P_U = ICF_P_U(eta,xa,P); 
      
        ckf = CKF(eta,xa,P); 

        for e=1:100 
            
            z = zt{e};
            zCount = zCountt{e}
            t_e = cputime;

            rgicf.prepData(z,zCount,H,Rinv);
            rgicf.consensus(K,E);
            rgicf.estimate(e,Phi,Q);     
            time_rgicf(nn) = time_rgicf(nn) + cputime - t_e;
            t_e = cputime;
            ggicf.prepData(z,zCount,H,Rinv);
            ggicf.consensus(K,E);
            ggicf.estimate(e,Phi,Q);
            
            time_ggicf(nn) = time_ggicf(nn) + cputime - t_e;
            t_e = cputime;

            sggicf.prepData(z,zCount,H,Rinv);
            sggicf.consensus(K,E,threshold);
            sggicf.estimate(e,Phi,Q);
            
            time_sggicf(nn) = time_sggicf(nn) + cputime - t_e;t_e = cputime;
            icf.prepData(z,zCount,H,Rinv);
            icf.consensus(K,eps,E);
            icf.estimate(e,Phi,Q);
            time_icf(nn) = time_icf(nn) + cputime - t_e; t_e = cputime;
            icf_P_U.prepData(z,zCount,H,Rinv);
            icf_P_U.consensus(K,eps,E,Tr,e);
            icf_P_U.estimate(e,Phi,Q);
            
            time_icf_P_U(nn) = time_icf_P_U(nn) + cputime - t_e;
            ckf.estimate(e,Phi,Q,z,zCount,H,Rinv);

        end 
        [  ME_rgicf(nn),SDE_rgicf(nn),...
            ME_ggicf(nn),SDE_ggicf(nn),...
            ME_sggicf(nn),SDE_sggicf(nn),...
            ME_ckf(nn),SDE_ckf(nn),...
            ME_icf(nn),SDE_icf(nn),...
            ME_icf_P_U(nn),SDE_icf_P_U(nn)]...
            = computeStats(xa,rgicf,ggicf,sggicf,ckf,icf,icf_P_U);
  
    end
%end
%%
figure(1)
plot(xa(1,1:e),xa(2,1:e),'k-','LineWidth',2); hold on
plot(xa(1,1),xa(2,1),'k<','LineWidth',2); hold on
set(legend,'FontSize',12);
grid on
axis on
xlabel('x [m]')
ylabel('y [m]')
title('Target')
xlim([0 1100])
ylim([0 1100])



figure(2)
h_gt=plot([1:e]*0.05, (xa(1,1:e)),'k-','LineWidth',1); hold on
h_icf=plot([1:e]*0.05, (icf(i).x(1,1:e)),'b-.','LineWidth',1.5);
hold on
h_icf_p_u=plot([1:e]*0.05, (icf_P_U(i).x(1,1:e)),'r--','LineWidth',2);
hold on
legend('target','ICF_{mod}','ICF');
grid on
axis on
xlabel('e')
ylabel('x [m]')
set(legend,'FontSize',12);


figure(3)
h_gt=plot([1:e]*0.05, (xa(2,1:e)),'k-','LineWidth',1); hold on
h_icf=plot([1:e]*0.05, (icf(i).x(2,1:e)),'b-.','LineWidth',1.5);
hold on
h_icf_p_u=plot([1:e]*0.05, (icf_P_U(i).x(2,1:e)),'r--','LineWidth',2);
hold on
legend('y-axis','ICF_{mod}','ICF');
grid on
axis on
xlabel('e ')
ylabel('y [m]')
set(legend,'FontSize',12);


%
figure(4);
m1 = sqrt(4);
for i=1:4
    if mod(i,m1)==0
        pos(i,1)=500;
        pos(i,2)=(floor(i/m1)-1)/(m1-1)*500;
    else
        pos(i,1)=(mod(i,m1)-1)/(m1-1)*500;
        pos(i,2)=floor(i/m1)/(m1-1)*500;
    end
    plot(pos(i,1),pos(i,2),'^','Markersize',8); hold on
end
E=zeros(4,4);
for ii=1:4
    if ii<=4-1
        E(ii,ii+1)=1;
    end
    if ii>=2
        E(ii,ii-1)=1;
    end
end
for  ii=1:4
    for jj=ii:4
       if E(ii,jj)==1
         plot([pos(ii,1),pos(jj,1)],[pos(ii,2),pos(jj,2)],'k-.','Linewidth',1); hold on
       end
    end
end
title('sensors')
xlim([-100 600])
ylim([-100 600])