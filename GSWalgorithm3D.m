

loadings= load('octa_load.txt','-ascii');
loading=loadings(1:end-1);% 0,1 only
trig=loadings(end);% 0,1 only
posXpath=gpuArray(single(load('octa_posXpaths.txt')));%single real
posYpath=gpuArray(single(load('octa_posYpaths.txt')));%single real
posZpath=gpuArray(single(load('octa_posZpaths.txt')));%single real
loading=gpuArray(uint32(loading));
trig=uint32(trig);
assigned=uint32(load('octa_assign.txt'));%positive integer
weights=gpuArray(single(ones(length(loading),1)));
peakvalues=gpuArray(single(zeros(length(loading),1)));%single complex
[~, framenumber]=(size(posXpath));

%% just local variables.
SLM_Xnum=512; SLM_Ynum=512;
CCD_num=512;
focal_length = 0.004;
lambda = 820E-9;
SLM_pixel = 15E-6;
telescope = 2;
beam_radius = 4E-3;
posMetric=lambda*focal_length*telescope/SLM_pixel/CCD_num;
%0.0012;%% 

% SLM_Xnum=800; SLM_Ynum=800;
% focal_length = 0.004;
% lambda = 850E-9;
% SLM_pixel = 16E-6;
% telescope = 3;
% beam_radius = 4E-3;
% posMetric=1;


posXpath=gpuArray(single(posXpath*posMetric));
posYpath=gpuArray(single(posYpath*posMetric));
posZpath=gpuArray(single(posZpath*posMetric));
xx=[1:SLM_Xnum]-SLM_Xnum/2-1;
[X, Y]=(meshgrid(xx,xx));
X=gpuArray(single(X));%single real
Y=gpuArray(single(Y));%single real

suminit1=gpuArray(single(zeros(SLM_Xnum)));%single complex
randVals=gpuArray(single(2*pi*rand(length(loading),1)));
suminit1=v2w(exp(1i*randVals),posXpath(:,1),posYpath(:,1),posZpath(:,1));

phaseGSW1=gpuArray(single(zeros(SLM_Xnum)));%single real
phaseGSW1=angle(suminit1);%(-pi : pi) domain 
phaseGSW11=gpuArray(uint32(zeros(SLM_Xnum)));% integer
phaseGSW11=mod(uint32((phaseGSW1+2*pi)*255/2/pi),256);% 0-255 domain. Matlab automatically rounding-off them.
phaseGSW2=gpuArray(single(zeros(512)));%single complex
phaseGSW2=cast(phaseGSW11,'single')*2*pi/255;
%phaseamp=gpuArray(single(ones(512)));%single real

c=2*beam_radius/SLM_pixel/2.35;
gaussianform=exp(-(X.*X+Y.*Y)/2/c/c);
circle=1-floor((X.*X+Y.*Y)/SLM_Xnum/SLM_Ynum*4);
circle=circle>0;
maske=gpuArray(gaussianform.*circle);%single real
% numover=4;%integer
% [XX,YY]=meshgrid(1:512*numover,1:512*numover);
% XX=gpuArray(single(XX));%single real
% YY=gpuArray(single(YY));%single real
% gaussianform=exp(-((XX-512*numover/2).^2+(YY-512*numover/2).^2)/2/118^2);%single real
% circle=1-floor(((XX-512*numover/2).^2+(YY-512*numover/2).^2)/(512/2)^2);%single real
% circle=circle>0;
% maske=gpuArray(gaussianform.*circle);%single real

V_in=gpuArray(single(ones(length(loading),1)));
targetAmp_pointwise=gpuArray(single(ones(length(loading),1)));
W_in=gpuArray(single(zeros(SLM_Ynum,SLM_Xnum)));
%W_in=suminit1;
% phaseover=gpuArray(single(ones(512*numover)));%single complex
% targetamp=gpuArray(single(zeros(512*numover)));%single real
% phaseR=gpuArray(single(zeros(512*numover)));%single complex

GSWiter=5;
rf=0.7;
sf=1;
psf=gpuArray(single(exp(-((X).^2/rf+(Y).^2/rf).^sf)));%single real
psf2=fft2(psf);%single complex
psf2=abs(psf2)/max(max(abs(psf2)));%to single real
% psf2=gpuArray(single(ones(512)));

weighter=gpuArray(single(zeros(framenumber,1))); %%
weighter(1:3)=[0.7 0.4 0.2]*0.01; weighter(4:framenumber-5)=0.1*0.01; weighter(framenumber-4:framenumber)=[0.1 0.1 0.4 0.7 1]*0.01;
    
    
%% GSW operation
% while 1; % run indefinitely. 
    
while trig==0;
    loadings= load('octa_load.txt');
    loading=loadings(1:end-1);
    trig=loadings(end);
end;
posXpath=gpuArray(single(load('octa_posXpaths.txt')));
posYpath=gpuArray(single(load('octa_posYpaths.txt')));
posZpath=gpuArray(single(load('octa_posZpaths.txt')));
%posXpath=reshape(posXpath,[length(loading),length(posXpath)/length(loading)]);
%posYpath=reshape(posYpath,[length(loading),length(posYpath)/length(loading)]);
posXpath=gpuArray(single(posXpath*posMetric));
posYpath=gpuArray(single(posYpath*posMetric));
posZpath=gpuArray(single(posZpath*posMetric));
loading=gpuArray(uint32(loading));
trig=uint32(trig);
assigned=uint32(load('octa_assign.txt'));
weights=gpuArray(single(ones(length(loading),1)));
peakvalues=gpuArray(single(zeros(length(loading),1)));
[~, framenumber]=(size(posXpath));



for qqq=1:framenumber;
     
    for gswiter=1:GSWiter*1;%[qqq, gswiter]
        V_in=w2v(gaussianform.*exp(1i*phaseGSW2),posXpath(:,qqq),posYpath(:,qqq),posZpath(:,qqq));
        peakvalues=V_in;
        
        if ((sum(loading)==0) || (sum(loading)==length(loading)))
            weights=weights.*mean(abs(peakvalues))./abs(peakvalues)/mean(weights);
            targetAmp_pointwise=weights;
        else
            if qqq>3 && qqq<framenumber-2
                weights(loading==0)=weighter(qqq);
            else
                weights(loading==1)=weights(loading==1).*mean(abs(peakvalues(loading==1)))./abs(peakvalues(loading==1))/mean(weights(loading==1));
                weights(loading==0)=weights(loading==0).*mean(abs(peakvalues(loading==0)))./abs(peakvalues(loading==0))/mean(weights(loading==0));
                weights(loading==0)=weights(loading==0)*weighter(qqq);
            end
            targetAmp_pointwise(assigned)=weights(assigned);
        end
%         for q=1:length(loading);% This loop is the bottleneck, 80% of execution time.
%             qq=assigned(q);
%             targetAmp_pointwise(qq)=weights(qq);
%             %targetamp(round(posYpath(qq,qqq)/posMetric/(1.271*8/numover)+512*numover/2),round(posXpath(qq,qqq)/posMetric/(1.271*8/numover)+512*numover/2))=weights(qq);
%         end;
        %% GSW algorithm in Wikipedia
        W_in=v2w(targetAmp_pointwise.*exp(1i*angle(V_in)),posXpath(:,qqq),posYpath(:,qqq),posZpath(:,qqq));
        phaseGSW1=angle(W_in);
        phaseGSW11=mod(uint32((phaseGSW1+2*pi)*255/2/pi),256);
        phaseGSW2=cast(phaseGSW11,'single')*2*pi/255;
%         phaseR=(ifft2(fftshift(targetamp.*exp(1i*angle(phaseoverfft)))));
%         phaseGSW1=angle(phaseR(512*numover/2-255:512*numover/2+256,512*numover/2-255:512*numover/2+256));
%         phaseGSW11=mod(uint32((phaseGSW1+2*pi)*255/2/pi),256);
%         figure; image(abs(plane(maske.*exp(1i*angle(W_in)),posZpath(1,1))).^2,'CDataMapping','scaled');
%         figure; image(abs(plane(maske.*exp(1i*angle(W_in)),posZpath(2,1))).^2,'CDataMapping','scaled');
   end;
    %phaseGSW111(:,:,qqq)=phaseGSW11;
    %     for q=1:length(loading);
    %         targetamp(round(posYpath(q,qqq)/(1.271*8/numover)+512*numover/2),round(posXpath(q,qqq)/(1.271*8/numover)+512*numover/2))=gpuArray(single(0));
    %     end;
    targetAmp_pointwise=gpuArray(single(ones(length(loading),1))); % either upper 'for' loop or this line, whichever faster.

    %% To be added by the user.
    % Call SLM library
    %
    % end;
     
%     points
%     trap_1(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-500*posMetric));
%     trap_2(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),0*posMetric));
%     trap_3(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),500*posMetric));

%     cube,sphere
%     trap_1(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-50*posMetric));
%     trap_2(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-30*posMetric));
%     trap_3(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-10*posMetric));
%     trap_4(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),10*posMetric));
%     trap_5(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),30*posMetric));
%     trap_6(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),50*posMetric));

%   DNA
%     trap_1(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-1000*posMetric));
%     trap_2(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-800*posMetric));
%     trap_3(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-600*posMetric));
%     trap_4(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-400*posMetric));
%     trap_5(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-200*posMetric));
%     trap_6(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),0*posMetric));
%     trap_7(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),200*posMetric));
%     trap_8(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),400*posMetric));
%     trap_9(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),600*posMetric));
%     trap_10(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),800*posMetric));
%     trap_11(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),1000*posMetric));
    
%   atom_orbit
%     trap_1(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-500*posMetric));
%     trap_2(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-400*posMetric));
%     trap_3(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-300*posMetric));
%     trap_4(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-200*posMetric));
%     trap_5(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),-100*posMetric));
%     trap_6(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),0*posMetric));
%     trap_7(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),100*posMetric));
%     trap_8(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),200*posMetric));
%     trap_9(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),300*posMetric));
%     trap_10(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),400*posMetric));
%     trap_11(:,:,qqq)=gather(plane(maske.*exp(1i*angle(W_in)),500*posMetric));
    
end;
% points
% k=1; figure; imagesc(abs(trap_1(:,:,k)).^2); figure; imagesc(abs(trap_2(:,:,k)).^2); figure; imagesc(abs(trap_3(:,:,k)).^2);

% cube, sphere
% trap_i=cat(3,trap_1(:,:,1),trap_2(:,:,1),trap_3(:,:,1),trap_4(:,:,1),trap_5(:,:,1),trap_6(:,:,1));
% trap_f=cat(3,trap_1(:,:,20),trap_2(:,:,20),trap_3(:,:,20),trap_4(:,:,20),trap_5(:,:,20),trap_6(:,:,20));
% figure; contourslice(abs(trap_i).^2,[],[],1:6); view(3); pbaspect([1 1 1]);
% figure; contourslice(abs(trap_f).^2,[],[],1:6); view(3); pbaspect([1 1 1]);
%%k=1; figure; imagesc(abs(trap_1(:,:,k)).^2); figure; imagesc(abs(trap_2(:,:,k)).^2); figure; imagesc(abs(trap_3(:,:,k)).^2);figure; imagesc(abs(trap_4(:,:,k)).^2);figure; imagesc(abs(trap_5(:,:,k)).^2);figure; imagesc(abs(trap_6(:,:,k)).^2);
% for j=1
%     figure; imagesc(abs(trap_1(:,:,j)).^2); set(gca,'Ydir','normal'); axis([150 350 250 450]);
%     title(['z=-62.3um ' int2str(j) 'th frame']);
% %     saveas(gca,['C:\Users\hknho\Desktop\3dgs\cube_simulation\cube' int2str(1) '-' int2str(j) '.png']);
% %     close;
%     
%     figure; imagesc(abs(trap_2(:,:,j)).^2); set(gca,'Ydir','normal'); axis([150 350 250 450]);
%     title(['z=-37.4um ' int2str(j) 'th frame']);
% %     saveas(gca,['C:\Users\hknho\Desktop\3dgs\cube_simulation\cube' int2str(2) '-' int2str(j) '.png']);
% %     close;
%     
%     figure; imagesc(abs(trap_3(:,:,j)).^2); set(gca,'Ydir','normal'); axis([150 350 250 450]);
%     title(['z=-12.5um ' int2str(j) 'th frame']);
% %     saveas(gca,['C:\Users\hknho\Desktop\3dgs\cube_simulation\cube' int2str(3) '-' int2str(j) '.png']);
% %     close;
%     
%     figure; imagesc(abs(trap_4(:,:,j)).^2); set(gca,'Ydir','normal'); axis([150 350 250 450]);
%     title(['z=12.5um ' int2str(j) 'th frame']);
% %     saveas(gca,['C:\Users\hknho\Desktop\3dgs\cube_simulation\cube' int2str(4) '-' int2str(j) '.png']);
% %     close;
%     
%     figure; imagesc(abs(trap_5(:,:,j)).^2); set(gca,'Ydir','normal'); axis([150 350 250 450]);
%     title(['z=37.4um ' int2str(j) 'th frame']);
% %     saveas(gca,['C:\Users\hknho\Desktop\3dgs\cube_simulation\cube' int2str(5) '-' int2str(j) '.png']);
% %     close;
%     
%     figure; imagesc(abs(trap_6(:,:,j)).^2); set(gca,'Ydir','normal'); axis([150 350 250 450]);
%     title(['z=62.3um ' int2str(j) 'th frame']);
% %     saveas(gca,['C:\Users\hknho\Desktop\3dgs\cube_simulation\cube' int2str(6) '-' int2str(j) '.png']);
% %     close;
% end


% % sphere
% % trap_i=cat(3,trap_1(:,:,1),trap_2(:,:,1),trap_3(:,:,1),trap_4(:,:,1),trap_5(:,:,1),trap_6(:,:,1));
% % trap_f=cat(3,trap_1(:,:,20),trap_2(:,:,20),trap_3(:,:,20),trap_4(:,:,20),trap_5(:,:,20),trap_6(:,:,20));
% % figure; contourslice(abs(trap_i).^2,[],[],1:6); view(3); pbaspect([1 1 1]);
% % figure; contourslice(abs(trap_f).^2,[],[],1:6); view(3); pbaspect([1 1 1]);
% %%k=1; figure; imagesc(abs(trap_1(:,:,k)).^2); figure; imagesc(abs(trap_2(:,:,k)).^2); figure; imagesc(abs(trap_3(:,:,k)).^2);figure; imagesc(abs(trap_4(:,:,k)).^2);figure; imagesc(abs(trap_5(:,:,k)).^2);figure; imagesc(abs(trap_6(:,:,k)).^2);
% figure;

% DNA
% trap_i=cat(3,trap_1(:,:,1),trap_2(:,:,1),trap_3(:,:,1),trap_4(:,:,1),trap_5(:,:,1),trap_6(:,:,1),trap_7(:,:,1),trap_8(:,:,1),trap_9(:,:,1),trap_10(:,:,1),trap_11(:,:,1));
% trap_f=cat(3,trap_1(:,:,20),trap_2(:,:,20),trap_3(:,:,20),trap_4(:,:,20),trap_5(:,:,20),trap_6(:,:,20),trap_7(:,:,20),trap_8(:,:,20),trap_9(:,:,20),trap_10(:,:,20),trap_11(:,:,20));
% figure; contourslice(abs(trap_i).^2,[],[],1:10); view(3); pbaspect([1 1 5]);
% figure; contourslice(abs(trap_f).^2,[],[],1:10); view(3); pbaspect([1 1 5]);
% %k=1; figure; imagesc(abs(trap_1(:,:,k)).^2); figure; imagesc(abs(trap_2(:,:,k)).^2); figure; imagesc(abs(trap_3(:,:,k)).^2);figure; imagesc(abs(trap_4(:,:,k)).^2);figure; imagesc(abs(trap_5(:,:,k)).^2);figure; imagesc(abs(trap_6(:,:,k)).^2);figure; imagesc(abs(trap_7(:,:,k)).^2);figure; imagesc(abs(trap_8(:,:,k)).^2);figure; imagesc(abs(trap_9(:,:,k)).^2);figure; imagesc(abs(trap_10(:,:,k)).^2);figure; imagesc(abs(trap_11(:,:,k)).^2);

% atom_orbit
% trap_i=cat(3,trap_1(:,:,1),trap_2(:,:,1),trap_3(:,:,1),trap_4(:,:,1),trap_5(:,:,1),trap_6(:,:,1),trap_7(:,:,1),trap_8(:,:,1),trap_9(:,:,1),trap_10(:,:,1),trap_11(:,:,1));
% trap_f=cat(3,trap_1(:,:,20),trap_2(:,:,20),trap_3(:,:,20),trap_4(:,:,20),trap_5(:,:,20),trap_6(:,:,20),trap_7(:,:,20),trap_8(:,:,20),trap_9(:,:,20),trap_10(:,:,20),trap_11(:,:,20));
% figure; contourslice(abs(trap_i).^2,[],[],2:10); view(3); pbaspect([1 1 1]);
% figure; contourslice(abs(trap_f).^2,[],[],2:10); view(3); pbaspect([1 1 1]);
%k=1; figure; imagesc(abs(trap_1(:,:,k)).^2); figure; imagesc(abs(trap_2(:,:,k)).^2); figure; imagesc(abs(trap_3(:,:,k)).^2);figure; imagesc(abs(trap_4(:,:,k)).^2);figure; imagesc(abs(trap_5(:,:,k)).^2);figure; imagesc(abs(trap_6(:,:,k)).^2);figure; imagesc(abs(trap_7(:,:,k)).^2);figure; imagesc(abs(trap_8(:,:,k)).^2);figure; imagesc(abs(trap_9(:,:,k)).^2);figure; imagesc(abs(trap_10(:,:,k)).^2);figure; imagesc(abs(trap_11(:,:,k)).^2);

