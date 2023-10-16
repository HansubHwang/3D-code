function V=plane(w_in,kZ)
    %w_in=W_in; kZ=posZpath(1,1);
    SLM_Xnum = 512; SLM_Ynum = 512;
    SLM_pixel = 16e-6; telescope = 3;
    lambda = 850E-9; focal_length = 0.004;
    
    X = repmat(gpuArray(([1:SLM_Xnum]-SLM_Xnum/2))*SLM_pixel/telescope,[SLM_Ynum,1]);
    Y = repmat(gpuArray(([1:SLM_Ynum]'-SLM_Ynum/2))*SLM_pixel/telescope,[1,SLM_Xnum]);
    
    V = gpuArray(single(zeros(SLM_Ynum,SLM_Xnum)));
    
%     numover = 4;
%     V0 = gpuArray(single(zeros(SLM_Ynum*numover,SLM_Xnum*numover)));
%     V0(512*numover/2-255:512*numover/2+256,512*numover/2-255:512*numover/2+256) = exp(-1i/lambda/focal_length/focal_length*kZ*(X.*X+Y.*Y)).*w_in;
%      V = fftshift(fft2(circshift(V0,[1 1])))/SLM_Xnum/SLM_Ynum;
     
    V = ifft2(exp(-1i*pi/lambda/focal_length/focal_length*kZ*(X.*X+Y.*Y)).*fft2(fftshift(fft2(w_in(1:SLM_Ynum,1:SLM_Xnum)))));
    
    