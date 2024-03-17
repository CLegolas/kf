classdef CKF < handle
   properties
       x;
        x_;
        J;
        J_;
   end 
    methods
        function a = CKF(eta,xa,P)
            if (nargin > 0)
                a = CKF;

                x_sum = zeros(4,1);
                a.J_ = zeros(4,4);
                for i=1:4
                    Jtemp = eye(4,4)/P{i};
                    a.J_ = a.J_ + Jtemp;
                    x_sum = x_sum + Jtemp * (xa(:,1) + eta{i});
                end
                a.x_ = a.J_\x_sum;
                
                
                a.x = zeros(4,100);
                a.J = zeros(4,4);
            end
        end 
  function estimate(a,t,Phi,Q,z,zCount,H,Rinv)
            u = zeros(4,4);
            U = zeros(4,4,4);
            for i=1:4
                if zCount(i)>0
                    u(:,i) = H' * Rinv * z{i};
                    U(:,:,i) = H' * Rinv * H;
                else
                    u(:,i) = zeros(4,1);
                    U(:,:,i) = zeros(4,4);                   
                end
            end
            
            x_sum = a.J_* a.x_;
            a.J = a.J_;
            for i=1:4
                a.J =  a.J + U(:,:,i);
                x_sum  = x_sum + u(:,i);
            end
            a.x(:,t) = a.J \ x_sum;
            a.x_ = Phi * a.x(:,t);
            a.J_ = eye(4,4)/ ( Phi / a.J * Phi' + Q );
            
end
end 
    
end