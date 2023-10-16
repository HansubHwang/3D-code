function W = v2w(v_in,kX,kY,kZ)
% W(xi,yj)=sum[|Vm|*exp(i*(angVm+2*pi/lambda/f*(kx,m*xi+ky,m*yj)+pi/lamda/f^2*kz)),m*(xi^2+yj^2)), for {m targets}]
%     v_in=([1:5]+1i*[1:5]).';
%     kX=[-400:200:400]'; kY=[400 200 0 200 400]'; kZ=[-10 -10 0 10 10]';
    %%constant
%     tic
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
    
%     toc
%     tic
    amp = repmat(abs(v_in)',[SLM_Xnum*SLM_Ynum,1]);
    phase = repmat(angle(v_in)',[SLM_Xnum*SLM_Ynum,1]) + 2*pi/lambda/focal_length*(X*kX'+Y*kY') + pi/lambda/focal_length/focal_length*((X.*X+Y.*Y)*kZ');
    W1 = amp.*exp(1i*phase);%
%     W2 = zeros(SLM_Xnum*SLM_Ynum,1);
    W2 = sum(W1,2);
%     toc
%     tic
%     for k = 1:length(v_in)
%         W2 = W2+W1(:,k);
%     end
%     toc
%     tic
    W=reshape(gpuArray(single(W2)),[SLM_Xnum,SLM_Ynum]);
%     toc
    
%     W=gpuArray(single(zeros(SLM_Ynum,SLM_Xnum)));
%     for n=1:length(kZ)
%        W=W+ifft2(sqrt(lambda*focal_length)*ifft2(exp(1i*pi*lambda*kZ(n)*(v_X.v_X+v_Y.v_Y)).*fft2(
%     end

