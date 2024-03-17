classdef SGICF < handle
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

        function a = SGICF(eta,xa,P)
            if (nargin > 0)
                a(4,1) = SGICF;
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
        function consensus(a,K,E,threshold)
            for k=1:K
                
                for i=1:4
                    a(i).vtemp = zeros(4,1);
                    a(i).Vtemp = zeros(4,4);    
                    for ij=1:4
                        fim(ij) = (a(ij).v - a(i).v)'/(a(ij).V + a(i).V)*(a(ij).v - a(i).v);
                    end  
                    neighbors = find(E(i,:));
                    nneighbors = numel(neighbors);
                    
                    prob = rand(1,nneighbors);
                  id = find(prob>threshold);
                    if isempty(id)
                        id = ceil(nneighbors*rand);
                    end
                    neighbor_temp = neighbors(id);
                    
                    [~,mi] = max(fim(neighbor_temp));
                    e = neighbor_temp(mi);
                    
                    a(i).vtemp = (a(i).v + a(e).v)/2;
                   a(i).Vtemp = (a(i).V + a(e).V)/2;
                    
                    a(i).v = a(i).vtemp;
                   a(e).v = a(i).vtemp;

                    a(i).V = a(i).Vtemp;
                    a(e).V = a(i).Vtemp;
                end
                
%                 for i=1:4
%                     
%                     a(i).v = a(i).vtemp;
%                     a(i).V = a(i).Vtemp;
%                     
%                 end
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