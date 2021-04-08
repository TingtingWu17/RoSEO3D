function FI=FI_weigth_calculation(signal,background,basis_matrix,secM)

%muz=sqrt(1-mux^2-muy^2);
%[muxx,muyy,muzz,muxy,muxz,muyz] = Quickly_rotating_matrix_muxyz_to_M(mux,muy,muz,gamma);

M = secM.';
I = signal*basis_matrix*M+background;  %the forward modle

% for out of plane or 3D case
FI = zeros(6,6);
FI(1,2) = signal^2*sum(basis_matrix(:,1).*basis_matrix(:,2)./I);
FI(1,3) = signal^2*sum(basis_matrix(:,1).*basis_matrix(:,3)./I);
FI(1,4) = signal^2*sum(basis_matrix(:,1).*basis_matrix(:,4)./I);
FI(1,5) = signal^2*sum(basis_matrix(:,1).*basis_matrix(:,5)./I);
FI(1,6) = signal^2*sum(basis_matrix(:,1).*basis_matrix(:,6)./I);

FI(2,3) = signal^2*sum(basis_matrix(:,2).*basis_matrix(:,3)./I);
FI(2,4) = signal^2*sum(basis_matrix(:,2).*basis_matrix(:,4)./I);
FI(2,5) = signal^2*sum(basis_matrix(:,2).*basis_matrix(:,5)./I);
FI(2,6) = signal^2*sum(basis_matrix(:,2).*basis_matrix(:,6)./I);

FI(3,4) = signal^2*sum(basis_matrix(:,3).*basis_matrix(:,4)./I);
FI(3,5) = signal^2*sum(basis_matrix(:,3).*basis_matrix(:,5)./I);
FI(3,6) = signal^2*sum(basis_matrix(:,3).*basis_matrix(:,6)./I);

FI(4,5) = signal^2*sum(basis_matrix(:,4).*basis_matrix(:,5)./I);
FI(4,6) = signal^2*sum(basis_matrix(:,4).*basis_matrix(:,6)./I);

FI(5,6) = signal^2*sum(basis_matrix(:,5).*basis_matrix(:,6)./I);

FI = FI+FI.';

FI(1,1) = signal^2*sum(basis_matrix(:,1).^2./I);
FI(2,2) = signal^2*sum(basis_matrix(:,2).^2./I);
FI(3,3) = signal^2*sum(basis_matrix(:,3).^2./I);
FI(4,4) = signal^2*sum(basis_matrix(:,4).^2./I);
FI(5,5) = signal^2*sum(basis_matrix(:,5).^2./I);
FI(6,6) = signal^2*sum(basis_matrix(:,6).^2./I);




end


function [muxx,muyy,muzz,muxy,muxz,muyz] = Quickly_rotating_matrix_muxyz_to_M(mux,muy,muz,gamma)
muxx = gamma.*mux.^2+(1-gamma)/3;
muyy = gamma.*muy.^2+(1-gamma)/3;
muzz = gamma.*muz.^2+(1-gamma)/3;
muxy = gamma.*mux.*muy;
muxz = gamma.*mux.*muz;
muyz = gamma.*muz.*muy;

end
