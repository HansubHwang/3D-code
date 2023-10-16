 function V = w2v(w_in,kX,kY,kZ)
% V=1/(SLM_Xnum*SLM_Ynum)*sum[W(xj,yk)*exp(-i*(angVm+2*pi/lambda/f*(kx,m*xi+ky,m*yj)+pi/lamda/f^2*kz)),for {SLM plane (xj,yk)}]

%      w_in=repmat(gpuArray(single(([0:0.1:51.1]+1i*[112.2:-0.2:10]))),[512,1]);
%      kX=gpuArray([-400:200:400]'); kY=gpuArray([400 200 0 200 400]'); kZ=gpuArray([-10 -10 0 10 10]');
%     tic
    % constant
    SLM_Xnum = 512; SLM_Ynum = 512;
    SLM_pixel = 16e-6; telescope = 3;
    lambda = 850E-9; focal_length = 0.004;
    
    v_X = [1:SLM_Xnum]';
    v_X = (v_X-SLM_Xnum/2)*SLM_pixel/telescope;
    v_Y = [1:SLM_Xnum]';
    v_Y = (v_Y-SLM_Xnum/2)*SLM_pixel/telescope;
    
    Xtemp = repmat(v_X',[SLM_Ynum,1]);
    X = reshape(Xtemp,[numel(Xtemp),1]);
    Y = repmat(v_Y,[SLM_Xnum,1]);
    
    amp = repmat(reshape(abs(w_in),[numel(w_in),1]),[1,length(kX)]);
    phi = reshape(angle(w_in),[numel(w_in),1]);
    phase = repmat(phi,[1,length(kX)]) - 2*pi/lambda/focal_length*(X*kX'+Y*kY') - pi/lambda/focal_length/focal_length*((X.*X+Y.*Y)*kZ'); %2*pi/512*round(1/lambda/focal_length*(X*kX'+Y*kY')*512)
    V1 = (amp.*exp(1i*phase)).';
    %V2 = zeros(length(kX),1);
    V2 = sum(V1,2);
    V = gpuArray(single(V2/SLM_Xnum/SLM_Ynum));
%     toc
%     W2 = reshape(w_in,[numel(w_in),1]);
%     phase = repmat(W2,[1,length(kX)]) - 2*pi/lambda/focal_length*(X*kX'+Y*kY') - pi/lambda/focal_length/focal_length*((X.*X+Y.*Y)*kZ');
%     V1 = (exp(1i*phase)).';
%     V2 = zeros(length(kX),1);
%     for k = 1:SLM_Xnum*SLM_Ynum
%         V2 = V2+V1(:,k);
%     end
%     V = gpuArray(single(V2.*exp(1i*2*pi/lambda*(2*focal_length+kZ))*SLM_pixel*SLM_pixel/lambda/focal_length/1i));
