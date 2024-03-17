classdef  ICF < handle
    properties
        x;
        x_;
        J;
        J_;
        V;
      Vtemp;
        v;
        vtemp;
    end
    methods
        function a = ICF(eta,xa,P)
            if (nargin > 0)
                a(4,1) = ICF;
                
                for i=1:4
                    a(i).x_= xa(:,1) + eta{i};
                    a(i).J_ = eye(4,4)/P{i};
                    a(i).x = zeros(4,100);
                  a(i).J = zeros(4,4);
                    a(i).v = zeros(4,1);

                    a(i).V = zeros(4,4);
                   a(i).vtemp = zeros(4,1);
                    a(i).Vtemp = zeros(4,4);
                end
            end
        end 
        function prepData(a,z,zCount,H,Rinv)
            for i=1:4
                if zCount(i)>0
                    u = H' * Rinv * z{i};
                    U = H' * Rinv * H;
                else
                    u = zeros(4,1);
                    U = zeros(4,4);                   
                end
                
               a(i).v = 1/4 * a(i).J_ * a(i).x_ +u;
                a(i).V = 1/4 * a(i).J_ + U;
            end
        end

        function consensus(a,K,eps,E)
            for k=1:K
                for i=1:4
                    a(i).vtemp = zeros(4,1);
              
                    a(i).Vtemp = zeros(4,4);
%                    Delta = sum(E(i,:));
                    for i_=1:4
                        if E(i,i_)
%                            a(i).vtemp = a(i).vtemp + epsi_ij_MW * ( a(i_).v - a(i).v ) ;
%                            a(i).Vtemp = a(i).Vtemp + epsi_ij_MW * ( a(i_).V - a(i).V ) ;   
                           a(i).vtemp = a(i).vtemp + eps * ( a(i_).v - a(i).v ) ;
                          a(i).Vtemp = a(i).Vtemp + eps * ( a(i_).V - a(i).V ) ;   
                        end
                    end
                    a(i).vtemp =   a(i).v  + a(i).vtemp;
                    a(i).Vtemp =   a(i).V  + a(i).Vtemp;
%                     eps =1;
%                     a(i).vtemp =   (1-Delta*eps)*a(i).v  + eps * a(i).vtemp;
%                     a(i).Vtemp =   (1-Delta*eps)*a(i).V  + eps * a(i).Vtemp;
                    a(i).v = a(i).vtemp;
                    a(i).V = a(i).Vtemp;
                end
            end
        end        
        function estimate(a,e,Phi,Q)
            for i=1:4
                a(i).x(:,e) = a(i).V \ a(i).v;
               a(i).J = 4 * a(i).V;
           
               a(i).x_ = Phi * a(i).x(:,e);
                a(i).J_ = eye(4,4)/ ( Phi / a(i).J * Phi' + Q );
 end
        end   
    end
end