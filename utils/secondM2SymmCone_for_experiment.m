function [mux,muy,muz,rotMobil] = secondM2SymmCone_for_experiment(secM,signal,background,basis_matrix, varargin)

s=opt2struct(varargin);
%% define a least square problem structure

% get the initial estimates
%------------------------------------------------------------
% construct the M matrix
M=[secM(1),secM(4),secM(5);....
    secM(4),secM(2),secM(6);...
    secM(5),secM(6),secM(3)];
[V,D]=eig(M);


% initialization via SVD
x0(1)=real(V(1,3));
x0(2)=real(V(2,3));
x0(3)=real(V(3,3));
x0(4)=1.5*real(D(3,3))-.5;


FI=FI_weigth_calculation(signal,background,basis_matrix,secM);

% check if weights for second moments are provided
% if isfield(s,'weight')
%     w=s.weight;
    w = FI;
    sizeSecM=size(secM);
    
    %objective function
    %------------------------------------------------------------
    
    objFcn=@(z)sum((symmCone2secM(z)'-secM)*w*(symmCone2secM(z)-secM'));
    
% else
%     objFcn=@(z)sum((symmCone2secM(z)-secM').^2);
% end



%  constraints in full 3D coordinate

lb=[-ones(3,1);zeros(1,1)];
ub=[ones(3,1);ones(1,1)];
%ub=[ones(3,1);100000];

options=optimoptions('fmincon','Display','none');

sol=fmincon(objFcn,x0,[],[],[],[],lb,ub,@mycon,options);

mux=sol(1);
muy=sol(2);
muz=sol(3);
rotMobil=sol(4);


    function out=symmCone2secM(z)
        z=reshape(z,1,4);
        muz_t=z(:,3);
        muxx=z(:,4).*z(:,1).^2+(1-z(:,4))/3;
        muyy=z(:,4).*z(:,2).^2+(1-z(:,4))/3;
        muzz=z(:,4).*muz_t.^2+(1-z(:,4))/3;
        muxy=z(:,4).*z(:,1).*z(:,2);
        muxz=z(:,4).*z(:,1).*muz_t;
        muyz=z(:,4).*z(:,2).*muz_t;
        out=[muxx,muyy,muzz,muxy,muxz,muyz]';
    end

    function [c,ceq] = mycon(z)
        c=[];
        ceq=sum(z(:,1:3).^2,2)-1;
    end
end

