
%% points
% [posXi,posYi]=meshgrid(-500:500:500,-500:500:500,-500:500:500); posZi=repmat([1000 1000 1000 1000 1000 1000 -1000 -1000 -1000 -1000],[5,1]);
% [posXf,posYf]=meshgrid(-300:200:300,-300:200:300); posZf=repmat([1000 1000 -1000 -1000 -1000],[10,1]);
% 
% posXi=reshape(posXi,[numel(posXi),1]); posYi=reshape(posYi,[numel(posYi),1]); posZi=reshape(posZi,[numel(posZi),1]);
% posXf=reshape(posXf,[numel(posXf),1]); posYf=reshape(posYf,[numel(posYf),1]); posZf=reshape(posZf,[numel(posZf),1]);
% 
% 
% 
% targeter=1:length(posXf);

%% cube
% [posXi,posYi,posZi]=meshgrid(-50:20:50,-50:20:50,-50:20:50);
% 
% posXi=reshape(posXi,[numel(posXi),1]); posYi=reshape(posYi,[numel(posYi),1]); posZi=reshape(posZi,[numel(posZi),1]);
% posYi=posYi+100;
% 
% targeter=[44:47 50:53 56:59 62:65 80:83 86:89 92:95 98:101 116:119 122:125 128:131 134:137 152:155 158:161 164:167 170:173];
% i_seq=[1 6 36 31 1 181 186 216 211 181 186 6 36 216 211 31];
% t_seq=[44 47 65 62 44 152 155 173 170 152 155 47 65 173 170 62];
% fg=figure; scatter3(posXi,posYi,posZi,repmat(20,[length(posXi),1])); hold on; scatter3(posXi(targeter),posYi(targeter),posZi(targeter),'filled'); plot3(posXi(i_seq),posYi(i_seq),posZi(i_seq),'b--'); plot3(posXi(t_seq),posYi(t_seq),posZi(t_seq),'r'); hold off; pbaspect([1 1 1]);
% saveas(fg,'C:\Users\hknho\Desktop\3dgs\cube_simulation\cube00.png');
%% sphere
% [X,Y,Z]=meshgrid(-500:100:500,-500:100:500,-500:200:500); sphere_index=(X.^2+Y.^2+Z.^2<=500^2);
% posXi=X(sphere_index); posYi=Y(sphere_index); posZi=Z(sphere_index);
% 
% posXi=reshape(posXi,[numel(posXi),1]); posYi=reshape(posYi,[numel(posYi),1]); posZi=reshape(posZi,[numel(posZi),1]);
% 
% seq=1:length(posXi); target=(posXi.^2+posYi.^2+posZi.^2<350^2);
% targeter=seq(target~=0);

%% tube
% N=10;
% theta1=2*pi/(N/2):2*pi/(N/2):2*pi;
% theta2=2*pi/N:2*pi/N:2*pi;
% posXi=repmat([100*cos(theta1),200*cos(theta1),300*cos(theta2),400*cos(theta2),500*cos(theta2)]',[6,1]);
% posYi=repmat([100*sin(theta1),200*sin(theta1),300*sin(theta2),400*sin(theta2),500*sin(theta2)]',[6,1]);
% posZi=repmat([-500:200:500],[2*N/2+3*N,1]);
% 
% posXi=reshape(posXi,[numel(posXi),1]); posYi=reshape(posYi,[numel(posYi),1]); posZi=reshape(posZi,[numel(posZi),1]);
% 
% seq=1:length(posXi); target=(posXi.^2+posYi.^2<=400^2)&(posXi.^2+posYi.^2>=300^2);
% targeter=seq(target~=0);

%% DNA
% N=11; theta=[0:2*pi/N:2*pi*(N-1)/N]';
% R=[-300:30:300];%[300,200,160,120,80,-80,-120,-160,-200,-300]';
% posXi=[cos(theta)*R];
% posYi=[sin(theta)*R];
% posZi=repmat([-1000:2000/(N-1):1000],[1,length(R)]); 
% 
% posXi=reshape(posXi,[numel(posXi),1]); posYi=reshape(posYi,[numel(posYi),1]); posZi=reshape(posZi,[numel(posZi),1]);
% 
% targeter=[1:33 67:77 111:121 155:165 199:231];%[1:33 56:66 89:99 133:143 166:176 199:231];

%% atom orbit
% N1=1; phi=2*pi/N1:2*pi/N1:2*pi; theta=acos(repmat(-500,[1,N1])/500);
% N2=6; phi=cat(2,phi,2*pi/N2:2*pi/N2:2*pi); theta=cat(2,theta,acos(repmat(-400,[1,N2])/500));
% N3=12; phi=cat(2,phi,2*pi/N3:2*pi/N3:2*pi); theta=cat(2,theta,acos(repmat(-300,[1,N3])/500));
% N4=24; phi=cat(2,phi,2*pi/N4:2*pi/N4:2*pi); theta=cat(2,theta,acos(repmat(-200,[1,N4])/500));
% N5=42; phi=cat(2,phi,2*pi/N5:2*pi/N5:2*pi); theta=cat(2,theta,acos(repmat(-100,[1,N5])/500));
% N6=12; phi=cat(2,phi,2*pi/N6:2*pi/N6:2*pi); theta=cat(2,theta,acos(repmat(0,[1,N6])/500));
% N7=42; phi=cat(2,phi,2*pi/N7:2*pi/N7:2*pi); theta=cat(2,theta,acos(repmat(100,[1,N7])/500));
% N8=24; phi=cat(2,phi,2*pi/N8:2*pi/N8:2*pi); theta=cat(2,theta,acos(repmat(200,[1,N8])/500));
% N9=12; phi=cat(2,phi,2*pi/N9:2*pi/N9:2*pi); theta=cat(2,theta,acos(repmat(300,[1,N9])/500));
% N10=6; phi=cat(2,phi,2*pi/N10:2*pi/N10:2*pi); theta=cat(2,theta,acos(repmat(400,[1,N10])/500));
% N11=1; phi=cat(2,phi,2*pi/N11:2*pi/N11:2*pi); theta=cat(2,theta,acos(repmat(500,[1,N11])/500));


% posXi=cat(2,500*cos(phi).*sin(theta),[-20 -20 20 20 -20 -20 20 20 0]); 
% posYi=cat(2,500*sin(phi).*sin(theta),[-20 20 -20 20 -20 20 -20 20 0]); 
% posZi=cat(2,500*cos(theta),[-20 -20 -20 -20 20 20 20 20 0]); 
% 
% posXi=reshape(posXi,[numel(posXi),1]); posYi=reshape(posYi,[numel(posYi),1]); posZi=reshape(posZi,[numel(posZi),1]);
% 
% seq=1:length(posXi); target=(abs(posZi-posYi*sqrt(3))<50)|(abs(posZi+posYi*sqrt(3))<50)|(abs(posZi)<50);
% targeter=seq(target~=0);


%% cube 3x3
% [posXi,posYi,posZi]=meshgrid(-30:30:30,-30:30:30,-30:30:30);
% 
% posXi=reshape(posXi,[numel(posXi),1]); posYi=reshape(posYi,[numel(posYi),1]); posZi=reshape(posZi,[numel(posZi),1]);
% 
% targeter=1:27;

%% octahedron_2plane
% a=10;
% posXi=[0  -sqrt(1.5) sqrt(1.5) 0 -sqrt(1.5) sqrt(1.5)]'*a;
% posYi=[sqrt(2) -1/sqrt(2) -1/sqrt(2) -sqrt(2) 1/sqrt(2) 1/sqrt(2)]'*a;
% posZi=[1 1 1 -1 -1 -1]'*a;
% posYi=posYi+10*a;
% targeter=1:6;

%% octahedron_3plane
% a=10;
% posXi=[0 -sqrt(1.5) -sqrt(1.5) sqrt(1.5) sqrt(1.5) 0]'*a;
% posYi=[0 -sqrt(1.5) sqrt(1.5) sqrt(1.5) -sqrt(1.5) 0]'*a;
% posZi=[-sqrt(3) 0 0 0 0 sqrt(3)]'*a;
% posYi=posYi+3*a;
% targeter=1:6;

%% 
% a=10;
% posXi=[-1 1 0 0]'*a;
% posYi=[0 0 -1 1]'*a;
% posZi=[-1 -1 1 1]'*a;
% posYi=posYi+3*a;
% targeter=1:4;

%%
AAA=linspace(-50,40,6); 
AAA2=linspace(-50,40,7);
posXi=[AAA' ; -50; 0; 1; 0;1; -16*sqrt(3)/2+2; -3+16*sqrt(3)/2; AAA2'];
posYi=[-60; -60; -60; -60; -60; -60; -20; -35; -17; -26; -8; -27; -15; 20; 20; 20; 20; 20; 20; 20];
posZi=zeros(20,1);
posZi([10,11,13:20]) = 9.5;
targeter=1:12;
%%
posXt=posXi(targeter); posYt=posYi(targeter); posZt=posZi(targeter);

% loading=randi(2,[length(posXi),1])-1;
loading=zeros(length(posXi),1);
% loading=ones(length(posXi),1);
for k=1:length(loading)
    if(rand(1)>0.5)
        loading(k)=1;
    end
end
% loading([1 3 4 7 9 10 15 16])=1;


costmat=zeros(length(posXt),length(posXi));
for kq=1:length(posXt);
    costmat(kq,:)=((posXi-posXt(kq)).^2+(posYi-posYt(kq)).^2+(posZi-posZt(kq)).^2)./loading;
end
costmat(isnan(costmat))=Inf;
assigned=munkres(costmat)';

targetsorted=1:length(posXt);
loading2=loading; loading2(assigned==0)=1;
assigned(assigned==0)=targeter(assigned==0); unassigned=setxor(1:length(loading),assigned); %%targetsorted->targeter

framenumber=20;
posXpath=zeros(length(posXi),framenumber); posYpath=zeros(length(posYi),framenumber); posZpath=zeros(length(posZi),framenumber);
for kqq=1:framenumber;
    posXpath(:,kqq)=posXi; posXpath(assigned,kqq)=posXi(assigned)-(posXi(assigned)-posXt)/(framenumber-1)*(kqq-1);
    posYpath(:,kqq)=posYi; posYpath(assigned,kqq)=posYi(assigned)-(posYi(assigned)-posYt)/(framenumber-1)*(kqq-1);
    posZpath(:,kqq)=posZi; posZpath(assigned,kqq)=posZi(assigned)-(posZi(assigned)-posZt)/(framenumber-1)*(kqq-1);
%     if kqq>framenumber/2
%         posXpath(loading2(targetsorted)==0,kqq)=posXi(assigned(assigned>length(posYt)));
%         posYpath(loading2(targetsorted)==0,kqq)=posYi(assigned(assigned>length(posXt)));
%         posZpath(loading2(targetsorted)==0,kqq)=posZi(assigned(assigned>length(posXt)));
%     end
end
save('octa_posXpaths.txt','posXpath','-ascii');
save('octa_posYpaths.txt','posYpath','-ascii');
save('octa_posZpaths.txt','posZpath','-ascii');
save('octa_assign.txt','unassigned','assigned','-ascii');
trig=1;
save('octa_load.txt','loading','trig','-ascii');

% % real position to [-649:650]
% posMetric=4.90e-7;
% posXpath= load('4_SHE_posXpaths.txt','-ascii')/posMetric;
% posYpath= load('4_SHE_posYpaths.txt','-ascii')/posMetric;
% posZpath= load('4_SHE_posZpaths.txt','-ascii')/posMetric;
% save('C:\Users\hknho\Desktop\v3\examples\ModifiedGSWAlgorithm 180405 backup\4_posXpaths.txt','posXpath','-ascii');
% save('C:\Users\hknho\Desktop\v3\examples\ModifiedGSWAlgorithm 180405 backup\4_posYpaths.txt','posYpath','-ascii');
% save('C:\Users\hknho\Desktop\v3\examples\ModifiedGSWAlgorithm 180405 backup\4_posZpaths.txt','posZpath','-ascii');