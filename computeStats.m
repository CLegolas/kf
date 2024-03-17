 function [ME_rgicf,SDE_rgicf,ME_ggicf,SDE_ggicf,ME_sggicf,SDE_sggicf,ME_ckf,SDE_ckf,ME_icf,SDE_icf,ME_icf_P_U,SDE_icf_P_U]...
    = computeStats(xa,rgicf,ggicf,sggicf,ckf,icf,icf_P_U)

e_rgicf = zeros(4,100);
e_ggicf = zeros(4,100);

e_sggicf = zeros(4,100);
e_ckf = zeros(1,100);
e_icf     = zeros(4,100);
e_icf_P_U = zeros(4,100);


for ee=1:100
    e_ckf(1,ee) = norm(xa(1:2,ee)-ckf.x(1:2,ee));
    for  i=1:4
        e_rgicf(i,ee)  = norm(xa(1:2,ee) -  rgicf(i).x(1:2,ee));
  e_ggicf(i,ee)  = norm(xa(1:2,ee) -  ggicf(i).x(1:2,ee));

        e_sggicf(i,ee) = norm(xa(1:2,ee) - sggicf(i).x(1:2,ee));

       e_icf(i,ee)     = norm(xa(1:2,ee) -     icf(i).x(1:2,ee));
        e_icf_P_U(i,ee) = norm(xa(1:2,ee) - icf_P_U(i).x(1:2,ee));
        
    end
end
ME_rgicf  = sum( e_rgicf(:))/400;
ME_ggicf  = sum( e_ggicf(:))/400;

ME_sggicf = sum(e_sggicf(:))/400;
ME_ckf = sum(e_ckf(:))/100; 

ME_icf      =    sum( e_icf(:))/400;
ME_icf_P_U = sum(e_icf_P_U(:))/400;
esqr_rgicf  = e_rgicf.*e_rgicf;
esqr_ggicf  = e_ggicf.*e_ggicf;

esqr_sggicf = e_sggicf.*e_sggicf;

esqr_ckf = e_ckf.*e_ckf;
esqr_icf  =     e_icf.*e_icf;
esqr_icf_P_U =  e_icf_P_U.*e_icf_P_U;
SDE_rgicf  = sqrt( sum(esqr_rgicf (:)) /400);
SDE_ggicf = sqrt( sum(esqr_ggicf(:)) /400);

SDE_sggicf = sqrt( sum(esqr_sggicf(:)) /400);
SDE_ckf = sqrt( sum(esqr_ckf(:)) /100);

SDE_icf =     sqrt( sum(esqr_icf(:)) /400);
SDE_icf_P_U = sqrt( sum(esqr_icf_P_U(:)) /400);

end