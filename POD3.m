function [Us,Vs,Ws,mkp,D] = POD3(U,V,W,METHOD,perc,keep_modes)
% FOR METHOD:
% 1 = Energy cutoff (requires a percentage input, else can enter [] for
% the cutoff percentage
% 2 = ELF (DCT-Entropy line fit) cutoff (see Brindise and Vlachos, 2017)
% 3 = Input desired modes for reconstruction (requires modes to be
% specified)

% If method is not specified, default is set to 2 (ELF) method
if nargin < 4
    METHOD = 2;
end

pc_NX = size(U,1);
pc_NY = size(U,2);
pc_NZ = size(U,3);

Umean = mean(U,4);
Vmean = mean(V,4);
Wmean = mean(W,4);
Us = U;
Vs = V;
Ws = W;

for t=1:size(U,4)
    Us(:,:,:,t) = U(:,:,:,t)-Umean;
    Vs(:,:,:,t) = V(:,:,:,t)-Vmean;
    Ws(:,:,:,t) = W(:,:,:,t)-Wmean;
end

tic,fprintf('POD filtering velocity errors...')

[M1,M2,M3,D,A] = KLDs3v(Us,Vs,Ws);

figure(17)
subplot(2,1,1)
loglog(D./sum(D),'-o')
title('Eigenspectra')
ylabel('Eigenvalue'); xlabel('Mode')
ylim([10^-10 10^0])
subplot(2,1,2)
semilogx(cumsum(D)./sum(D),'-o')
title('Psuedo Energy')
ylabel('Percent Energy'); xlabel('Mode')
% % % 
% % % for t=1:min(size(U,3),40)
% % %     figure(18)
% % %     subplot(10,4,t)
% % %     plot(A(t,:))
% % %     ylabel(num2str(t))
% % % end

% Moving average smoothing of temporal coefficients
As = A;
As(:,1        ) = 1/5*(A(:,1        ) + A(:,2        ) + A(:,3        ) + A(:,  size(U,4)-1) + A(:,  size(U,4)-0) );
As(:,2        ) = 1/5*(A(:,1        ) + A(:,2        ) + A(:,3        ) + A(:,4        ) + A(:,  size(U,4)-0) );
As(:,3:size(U,4)-2) = 1/5*(A(:,1:size(U,4)-4) + A(:,2:size(U,4)-3) + A(:,3:size(U,4)-2) + A(:,4:size(U,4)-1) + A(:,5:size(U,4)-0) );
As(:,  size(U,4)-1) = 1/5*(A(:,1        ) + A(:,  size(U,4)-3) + A(:,  size(U,4)-2) + A(:,  size(U,4)-1) + A(:,  size(U,4)-0) );
As(:,  size(U,4)  ) = 1/5*(A(:,1        ) + A(:,2        ) + A(:,  size(U,4)-2) + A(:,  size(U,4)-1) + A(:,  size(U,4)-0) );


for t=1:min(size(U,4),40)
    figure(18)
    subplot(10,4,t)
    plot(A(t,:),'b')
    hold on
    plot(As(t,:),'r')
    %imagesc(M1(:,:,t))
    %axis off
    title(['Mode ',num2str(t)]) 
    %ylabel(num2str(t))
    hold off
end

%%% SELECTION OF MODES TO USE IN RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use energy to determine modes to use
if METHOD == 1
    NN = find(cumsum(D)/sum(D)>perc,1,'first');
    %Remove first mode1
    %Dno1 = D(2:end);
    %NN = find(cumsum(Dno1)/sum(Dno1)>perc,1,'first');
    NN = max(NN,1);
    mkp = 1:NN;
end

%Use ELF method to determine modes to keep
if METHOD == 2
    for q = 1:1:length(D)
        m = M1(:,:,:,q); % U eigenmodes
        n = M2(:,:,:,q); % V eigenmodes
        p = M3(:,:,:,q); % W eigenmodes
        
        dctm = fftn(m); % Take FFT of U, V, and W eigenmodes in order to
        dctn = fftn(n); % compress eigenmodes to uniform coefficients.
        dctp = fftn(p); % For 3D, cannot use DCT
        dctm2 = abs(dctm).^2; % Take square of abs magnitude for entropy 
        dctn2 = abs(dctn).^2; % calculation.
        dctp2 = abs(dctp).^2;
        
        dctma2 = dctm2/sum(dctm2(:)); % Normalize DCT magnitudes so sum of 
        dctna2 = dctn2/sum(dctn2(:)); % field is 1 (for entropy calc)
        dctpa2 = dctp2/sum(dctp2(:));
        
        dctma2( dctma2 == 0 ) = 1; % Set all 0 values to 1 so they becomes
        dctna2( dctna2 == 0 ) = 1; % 0 in the entropy calculation sin log(1) = 0
        dctpa2( dctpa2 == 0 ) = 1;
        
        entm2(q) = sum(sum(sum(-log2(dctma2).*dctma2))); % Compute the shannon entropy
        entn2(q) = sum(sum(sum(-log2(dctna2).*dctna2))); % of the eigenmodes
        entp2(q) = sum(sum(sum(-log2(dctpa2).*dctpa2))); % of the eigenmodes
    end
    
    entmin = min([entm2;entn2;entp2],[],1); % Take the minimum between the U, V, and W eigenmodes
    [srtmin,indmin] = sort(entmin,'ascend'); % Sort by minimum entropy so modes out
                                             % of order now

    [bestind,lineinfo] = twolinefit_v2(srtmin); % Use two-line fit to determine the
    kpinds = indmin(1:bestind);                 % optimal cutoff mode
    
    figure(19); plot(srtmin,'k') % Plot sorted entropy
    hold on; plot([bestind,bestind],[0,max(entmin)],'r') % Plot location of cutoff
    save('entropy_line_fit.mat','srtmin')
    mkp = sort(kpinds,'ascend'); % Save the modes to keep based on optimal cutoff mode
end

% Use desired input set of modes for reconstruction
if METHOD == 3
    mkp = sort(keep_modes,'ascend');
end


%%% RECONSTRUCTION OF VELOCITY FIELDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('reconstructing %u modes',length(mkp))
fprintf(' and %.6g percent energy\n', sum(D(mkp))/sum(D)*100 )

%use NN timesteps with smoothed coefficients
% keyboard
Us = repmat(Umean,[1,1,1,size(U,4)])+reshape(reshape(M1(:,:,:,mkp),pc_NX*pc_NY*pc_NZ,length(mkp))*As(mkp,:),pc_NX,pc_NY,pc_NZ,size(U,4));
Vs = repmat(Vmean,[1,1,1,size(U,4)])+reshape(reshape(M2(:,:,:,mkp),pc_NX*pc_NY*pc_NZ,length(mkp))*As(mkp,:),pc_NX,pc_NY,pc_NZ,size(U,4));
Ws = repmat(Wmean,[1,1,1,size(U,4)])+reshape(reshape(M3(:,:,:,mkp),pc_NX*pc_NY*pc_NZ,length(mkp))*As(mkp,:),pc_NX,pc_NY,pc_NZ,size(U,4));


fprintf('%g\n',toc)

Us = Us;
Vs = Vs;
Ws = Ws;
