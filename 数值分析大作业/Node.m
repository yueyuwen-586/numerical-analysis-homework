classdef Node<handle
   properties
      isleaf
      a
      c
      invf
      invpsil
      invpsir
      alphal
      alphar
      betal
      betar
      deltal
      deltar
      miul
      miur
      sigma
      jl
      jr
   end
   properties (SetAccess = private)
      Next1 = Node.empty
      Next2 = Node.empty
      Prev = Node.empty
   end
   methods
       function f=LeafList(node)
           if node.isleaf==true
              f=node;
           else
              f=[LeafList(node.Next1) LeafList(node.Next2)];
           end
       end
       function solvesigma(node)
           node.sigma=node.miul*node.invpsil+node.miur*node.invpsir+node.invf;
       end
       function devide(node,c,K,f,g1,g2,psi1,psi2)
           if node.isleaf==false
               if node.Next1.isleaf==true && node.Next2.isleaf==true && monitor(node.Next1)+monitor(node.Next2)<c/(2^K)
                       cutNode(node);
               else
                   devide(node.Next1,c,K,f,g1,g2,psi1,psi2);
                   devide(node.Next2,c,K,f,g1,g2,psi1,psi2);
               end
           else 
               if monitor(node)>c
                   grownode(node,K,f,g1,g2,psi1,psi2);
               end
           end
       end
       function f=monitormax(node)
           if node.isleaf==true
               f=monitor(node);
           else
               f=max(monitormax(node.Next1),monitormax(node.Next2));
           end
       end
       function f=monitor(node)
           K=max(size(node.sigma));
           f=abs(node.sigma(K-1))+abs(node.sigma(K)-node.sigma(K-2));
       end
       function upsweepnumber(node,node1,node2)
           adl=node1.alphal;
           adr=node1.alphar;
           ael=node2.alphal;
           aer=node2.alphar;
           bdl=node1.betal;
           bdr=node1.betar;
           bel=node2.betal;
           ber=node2.betar;
           ddl=node1.deltal;
           ddr=node1.deltar;
           del=node2.deltal;
           der=node2.deltar;
           k=1-aer*bdl;
           node.alphal=(1-ael)*(adl-bdl*aer)/k+ael;
           node.alphar=aer*(1-bdr)*(1-adl)/k+adr;
           node.betal=bdl*(1-ber)*(1-ael)/k+bel;
           node.betar=(1-bdr)*(ber-bdl*aer)/k+bdr;
           node.deltal=(1-ael)*ddl/k+del+(ael-1)*bdl*der/k;
           node.deltar=(1-bdr)*der/k+ddr+(bdr-1)*aer*ddl/k;
       end
       function downsweepnumber(node,node1,node2)
           adl=node1.alphal;
           adr=node1.alphar;
           ael=node2.alphal;
           aer=node2.alphar;
           bdl=node1.betal;
           bdr=node1.betar;
           bel=node2.betal;
           ber=node2.betar;
           ddl=node1.deltal;
           ddr=node1.deltar;
           del=node2.deltal;
           der=node2.deltar;
           node1.miul=node.miul;
           node2.miur=node.miur;
           v=[1,aer;bdl,1]\[node.miur*(1-ber)-der;node.miul*(1-adl)-ddl];
           node1.miur=v(1);
           node2.miul=v(2);
       end
       function grownode(node,K,f,g1,g2,psi1,psi2)
           node.isleaf=false;
           b1=node.a;
           b2=node.c;
           b3=node.a+node.c;
           b3=b3/2;
           node.Next1=Node(b1,b3);
           node.Next2=Node(b3,b2);
           solveinv(node.Next1,K,f,g1,g2,psi1,psi2);
           evaluatenumber(node.Next1,g1,g2);
           solveinv(node.Next2,K,f,g1,g2,psi1,psi2);
           evaluatenumber(node.Next2,g1,g2);
       end
       function upsweep(node)
           if node.isleaf==false
               upsweep(node.Next1);
               upsweep(node.Next2);
               upsweepnumber(node,node.Next1,node.Next2);
           end
       end
       function downsweep(node)
           if node.isleaf==false
               downsweepnumber(node,node.Next1,node.Next2);
               downsweep(node.Next1);
               downsweep(node.Next2);
           end
       end
       function solveinv(node,K,f,g1,g2,psi1,psi2)
           node.invf=Pinverse(f,K,node.a,node.c,g1,g2,psi1,psi2);
           node.invpsil=Pinverse(psi1,K,node.a,node.c,g1,g2,psi1,psi2);
           node.invpsir=Pinverse(psi2,K,node.a,node.c,g1,g2,psi1,psi2);
       end
       function evaluatenumber(node,g1,g2)
           b1=node.a;
           b2=node.c;
           K=max(size(node.invf));
           for j=1:K
           h1=@(t) g1(t).*cos((j-1)*acos(2*(t-b1)/(b2-b1)-1));
           h2=@(t) g2(t).*cos((j-1)*acos(2*(t-b1)/(b2-b1)-1));
           a1=integral(h1,b1,b2);
           a2=integral(h2,b1,b2);
           %t=cheb(b1,b2,K);
           %v1=h1(t);
           %v2=h2(t);
           %w1=chebint(v1,b1,b2,K);
           %w2=chebint(v2,b1,b2,K);
           %a1=chebintegral(K,w1,b1,b2);
           %a2=chebintegral(K,w2,b1,b2);
           node.alphal=node.alphal+node.invpsil(j)*a1;
           node.alphar=node.alphar+node.invpsil(j)*a2;
           node.betal=node.betal+node.invpsir(j)*a1;
           node.betar=node.betar+node.invpsir(j)*a2;
           node.deltal=node.deltal+node.invf(j)*a1;
           node.deltar=node.deltar+node.invf(j)*a2;
           end
       end
      function node = Node(a,c)
         if (nargin > 0)
            node.isleaf = true;
            node.a=a;
            node.c=c;
            node.alphal=0;
            node.alphar=0;
            node.betal=0;
            node.betar=0;
            node.deltal=0;
            node.deltar=0;
            node.miul=0;
            node.miur=0;
         end
      end
      function insertBefore(newNode, nodeAfter1,nodeAfter2)
          newNode.Next1 = nodeAfter1;
          newNode.Next2 = nodeAfter2;
          nodeAfter1.Prev = newNode;
          nodeAfter2.Prev = newNode;
          newNode.isleaf=false;
      end
      function cutNode(node)
         delete(node.Next1);
         delete(node.Next2);
         node.Next1 =   Node.empty;
         node.Next2 =   Node.empty;
         node.isleaf=true;
      end

   end
    methods (Access = private)
      function delete(~)
      end
    end
end