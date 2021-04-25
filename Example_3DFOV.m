%%

% This 3D fabric orientation visualiation technique is designed for  
% various 3D dataset, here, the CT dataset is implemented as a example.
% For more detailed information, please refer to this journal paper:
% Jiang, S., et al., Three-dimensional fabric orientation visualisation
% technique for distributed fracture using X-ray computed tomography. 
% If you use this technique for academic purposes, please cite the abovementioned paper.


% the patch size is determined by patchw, which is 2*patchw. The patch size
% can also be the odd number. The overlap between patches is determined by the
% step. 

% For the larger patchw value, you should better increase the N number as
% well.


%%
% Determine the patch parameters
patchw=16; % patch size
step=16 ; % patch step

%% 
% Input the 3D dataset
for i=380:1:500 
       D(:,:,i-379)=imread(strcat('l_0',num2str(i),'.bmp'));
end


%%
% Divide the 3D dataset into overlapping patches
[m,n,p]=size(D);
gridx=patchw:step:m-patchw;
gridy=patchw:step:n-patchw; 
gridz=patchw:step:p-patchw;
D=double(D);


%%
% Compute a radial hamming window  
w=NaN.*zeros(patchw*2,patchw*2,patchw*2);
for k=1:patchw*2  
    for j=1:patchw*2  
        for i=1:patchw*2 
            dst=sqrt((i-0.5-patchw)^2+(j-0.5-patchw)^2+(k-0.5-patchw)^2) ; 
            w(i,j,k)=0.5*(cos(2*pi*dst/(patchw*2)))+0.5 ; 
            if (dst>patchw) 
                w(i,j,k)=0 ; 
            end 
        end
    end 
end 


%% 
% Compute the structure tensor mask_Q (MonteCarlo)
 
n_nbmaskQ=zeros(patchw*2,patchw*2,patchw*2) ;  
n_maskQ=zeros(patchw*2,patchw*2,patchw*2,3,3) ; 

for j=1:10000 % Number of iteration in Monte-Carlo simulation
N=10000 ; % Number of particle throw per iteration

r=(rand(N,3)-0.5)*2*patchw ;  
K(:,1)=r(:,1)./(sqrt(r(:,1).^2+r(:,2).^2+r(:,3).^2)) ;   
K(:,2)=r(:,2)./(sqrt(r(:,1).^2+r(:,2).^2+r(:,3).^2)) ; 
K(:,3)=r(:,3)./(sqrt(r(:,1).^2+r(:,2).^2+r(:,3).^2)) ;   

x=floor(r+patchw)+1 ;  
   for i=1:N
   n_nbmaskQ(x(i,1),x(i,2),x(i,3))=n_nbmaskQ(x(i,1),x(i,2),x(i,3))+1 ; % Count the number of particles of each voxel
   
   
   n_maskQ(x(i,1),x(i,2),x(i,3),1,1)=n_maskQ(x(i,1),x(i,2),x(i,3),1,1)+K(i,1)*K(i,1) ; 
   n_maskQ(x(i,1),x(i,2),x(i,3),1,2)=n_maskQ(x(i,1),x(i,2),x(i,3),1,2)+K(i,1)*K(i,2) ; 
   n_maskQ(x(i,1),x(i,2),x(i,3),1,3)=n_maskQ(x(i,1),x(i,2),x(i,3),1,3)+K(i,1)*K(i,3) ; 
   
   n_maskQ(x(i,1),x(i,2),x(i,3),2,1)=n_maskQ(x(i,1),x(i,2),x(i,3),2,1)+K(i,2)*K(i,1) ; 
   n_maskQ(x(i,1),x(i,2),x(i,3),2,2)=n_maskQ(x(i,1),x(i,2),x(i,3),2,2)+K(i,2)*K(i,2) ;
   n_maskQ(x(i,1),x(i,2),x(i,3),2,3)=n_maskQ(x(i,1),x(i,2),x(i,3),2,3)+K(i,2)*K(i,3) ; 
   
   n_maskQ(x(i,1),x(i,2),x(i,3),3,1)=n_maskQ(x(i,1),x(i,2),x(i,3),3,1)+K(i,3)*K(i,1) ; 
   n_maskQ(x(i,1),x(i,2),x(i,3),3,2)=n_maskQ(x(i,1),x(i,2),x(i,3),3,2)+K(i,3)*K(i,2) ;
   n_maskQ(x(i,1),x(i,2),x(i,3),3,3)=n_maskQ(x(i,1),x(i,2),x(i,3),3,3)+K(i,3)*K(i,3) ; 
   end

end

% Final scaling for the n_maskQ coefficients.    

% range (-1,1) scaling the possibility
n_maskQ(:,:,:,1,1)=n_maskQ(:,:,:,1,1)./n_nbmaskQ ; 
n_maskQ(:,:,:,1,2)=n_maskQ(:,:,:,1,2)./n_nbmaskQ ; 
n_maskQ(:,:,:,1,3)=n_maskQ(:,:,:,1,3)./n_nbmaskQ ; 
n_maskQ(:,:,:,2,1)=n_maskQ(:,:,:,2,1)./n_nbmaskQ ; 
n_maskQ(:,:,:,2,2)=n_maskQ(:,:,:,2,2)./n_nbmaskQ ; 
n_maskQ(:,:,:,2,3)=n_maskQ(:,:,:,2,3)./n_nbmaskQ ;
n_maskQ(:,:,:,3,1)=n_maskQ(:,:,:,3,1)./n_nbmaskQ ; 
n_maskQ(:,:,:,3,2)=n_maskQ(:,:,:,3,2)./n_nbmaskQ ; 
n_maskQ(:,:,:,3,3)=n_maskQ(:,:,:,3,3)./n_nbmaskQ ;


%% 
% calculate all the tensors
for k=1:length(gridz)  % Loop over the grids means each patch in the dataset
   for j=1:length(gridy) 
      for i=1:length(gridx) 

          
            % find out the centroid point of each patch
            % define a grid where the orientation matrix will be computed
            % Extract the single patch
            patch=D(gridx(i)-patchw+1:gridx(i)+patchw,gridy(j)-patchw+1:gridy(j)+patchw,gridz(k)-patchw+1:gridz(k)+patchw) ; 
            
            % pro-process of the patch data
            % Intensity scaling
            if std(patch(:))~=0  
                patch=(patch-mean(mean(mean(patch))))/std(patch(:)); 
            else
                patch=(patch-mean(mean(mean(patch)))) ;
            end 
           
                       
            S=fftshift(abs(fftn(patch.*w).^2)) ; % Fourier transform to get power spectrum (Fig.1)

            % The following compute the eq. 3 to obtain the gradient tensor
            % Just to handle properly empty patches
            if (sum(sum(sum(S)))==0) 
                Q(:,:,:,1,1)=S ; Q(:,:,:,1,2)=S ; Q(:,:,:,1,3)=S ;
                Q(:,:,:,2,1)=S ; Q(:,:,:,2,2)=S ; Q(:,:,:,2,3)=S ;
                Q(:,:,:,3,1)=S ; Q(:,:,:,3,2)=S ; Q(:,:,:,3,3)=S ;
            else
                % Actual usefull part  
                Q(:,:,:,1,1)=n_maskQ(:,:,:,1,1).*S./sum(sum(sum(S))) ; 
                Q(:,:,:,1,2)=n_maskQ(:,:,:,1,2).*S./sum(sum(sum(S))) ; 
                Q(:,:,:,1,3)=n_maskQ(:,:,:,1,3).*S./sum(sum(sum(S))) ; 
                
                Q(:,:,:,2,1)=n_maskQ(:,:,:,2,1).*S./sum(sum(sum(S))) ; 
                Q(:,:,:,2,2)=n_maskQ(:,:,:,2,2).*S./sum(sum(sum(S))) ; 
                Q(:,:,:,2,3)=n_maskQ(:,:,:,2,3).*S./sum(sum(sum(S))) ; 
                
                Q(:,:,:,3,1)=n_maskQ(:,:,:,3,1).*S./sum(sum(sum(S))) ; 
                Q(:,:,:,3,2)=n_maskQ(:,:,:,3,2).*S./sum(sum(sum(S))) ; 
                Q(:,:,:,3,3)=n_maskQ(:,:,:,3,3).*S./sum(sum(sum(S))) ; 
            end
                                  
            % Q2: 5D matrix containing the whole results. First 2 indices: structure tensor, indices 3 & 4 & 5: patch location,            
            
            Q2(:,:,i,j,k)=permute(sum(sum(sum(Q,1),2),3),[4,5,1,2,3]) ; 
            
            % Compute the anisotropic part with respect to eq.8
            T(:,:,i,j,k)=Q2(:,:,i,j,k);
            T(1,1,i,j,k)=Q2(1,1,i,j,k)-1/3 ;
            T(2,2,i,j,k)=Q2(2,2,i,j,k)-1/3 ; 
            T(3,3,i,j,k)=Q2(3,3,i,j,k)-1/3 ;
            T(:,:,i,j,k)=sqrt(3).*T(:,:,i,j,k)./sqrt(2) ;
            T_order(i,j,k)=(sum(sum(T(:,:,i,j,k).^2))).^(0.5);
      end
   end
end



%% 

% Plot the 3D ellipsoial visualisation of the raw dataset acorrding to the gradient
% tensors

figure ('Name','3D ellipsoid visualiation')


% Choose the appropriate colormap
C1=colormap('cool');
[m,n]=size(C1);

% Determine the range of the ordering parameters  in order to beter scale
% the ellipsoids

T_min=min(T_order(:));
T_max=max(T_order(:));
Len=T_max-T_min;


% Plot 3D images
for k=1 :length(gridz)  
   for j=1:length(gridy)
      for i=1:length(gridx)
%            Rescaling the colormap according to T_nor
           T_nor(i,j,k)=(T_order(i,j,k)-T_min)./Len;
           [newx(:,:,i,j,k), newy(:,:,i,j,k), newz(:,:,i,j,k), C]=Ellipsoidvisualisation(Q2(:,:,i,j,k),T_order(i,j,k),T_nor(i,j,k),i,j,k);         
           surf(newx(:,:,i,j,k), newy(:,:,i,j,k), newz(:,:,i,j,k),C,'FaceAlpha',T_nor(i,j,k)) ;                         
           shading interp;
           hold on;
      end
   end
end

caxis([0 max(T_order(:))])
axis equal  
view(45,45)
xlabel('X'),ylabel('Y'),zlabel('Z')  
h=colorbar;
set(h, 'YTick', [0 0.2 0.4 0.6 0.8 1]) 
set(h,'YTickLabel',{'0','0.2','0.4','0.6','0.8','1'})
set(h,'FontName','Times New Roman','FontSize',24)





