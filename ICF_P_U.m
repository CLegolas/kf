classdef ICF_P_U < handle
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
        function a = ICF_P_U(eta,xa,P)
            if (nargin > 0)
                a(4,1) = ICF_P_U;
                
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
        function consensus(a,K,eps,E,Tr,e)
           A=E;  
           D=zeros(4);
           for ii=1:4
               D(ii,ii)=sum( E(ii,:) );
           end
           L = D - A;
           eig_val = eig(L);                   
           eig_val = sort(eig_val,'ascend');   
%            eig_diff_val(1) = eig_val(1) ;
%            k=1 ;
%            for ii=1:4
%                if eig_val(ii) ~= eig_diff_val(end)
%                  k=k+1 ;
%                  eig_diff_val(k) = eig_val(ii) ;
%                end
%            end
           eig_diff_val = eig_val ;
           DD = length(eig_diff_val) -1;
           N_max = max( sum(A,2) ) ;
           for k=1:K
%                NN =  ceil(k/2) ;
%                if NN<=DD
%                    beta = 1 / eig_diff_val(DD+2-NN);
%                    ESI  = eye(4)-beta*L ;
%                 else
%                    ESI  = eye(4)-eps *  E;
%                end
               
               NN = mod(k,DD);
               if NN==0
                  NN=DD ;
               end
               if k<=DD*2
                   beta = 1 / eig_diff_val(DD+2-NN);
                   ESI  = eye(4)-beta*L ;
               else
                   ESI  = eye(4)-eps*E;
               end
               
                for i=1:4
                    a(i).vtemp = zeros(4,1);
                    a(i).Vtemp = zeros(4,4);
                         for i_=1:4
                         if E(i,i_)                         
                           Tri        = Tr(:,:,mod(k,2)+1) ;
                           
                           a(i).vtemp = a(i).vtemp + ESI(i,i_) * Tri * ( a(i_).v - a(i).v ) ;
                           a(i).Vtemp = a(i).Vtemp + ESI(i,i_) * Tri * ( a(i_).V - a(i).V ) ;   
                           
                           
                        end
                    end
                    a(i).vtemp =   a(i).v  + a(i).vtemp;
                    a(i).Vtemp =   a(i).V  + a(i).Vtemp;
                    
                end

%                     for i_=1:4
%                         if E(i,i_) 
% %                          a(i).vtemp =Tri* a(i_).v + ( eye(4,4) - Tri ) * a(i).v ;
% %                          a(i).Vtemp =Tri* a(i_).V + ( eye(4,4) - Tri ) * a(i).V ;                    
%                         end
%                     end
%                     a(i).vtemp =   a(i).v  + a(i).vtemp;
%                     a(i).Vtemp =   a(i).V  + a(i).Vtemp;
%                     
%                 end
                
                  for i=1:4
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