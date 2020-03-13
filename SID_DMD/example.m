% Written by Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (victor.zavala@wisc.edu)

rng(1)

close all
clear
clc

% create data 
N = 300;
l = 100;
A = [cos(pi/30) -sin(pi/30);sin(pi/30) cos(pi/30)];
snapind = [1,6,11,16,21,26,31,36,41,46,51,56,61];
X = zeros(2,N);
X(:,1) = [.9999999;0];
for i=2:N
    X(:,i) = A*X(:,i-1)+randn(2,1)*0;
end
wiscraw = [zeros(233,8) imread('wisc.png','png')==0 zeros(233,8)];
wisc = zeros(l);
wisc(l/5+1:l*4/5,l/5+1:l*4/5) =double(imresize(wiscraw,[l*3/5,l*3/5]));
C = [wisc(:) zeros(l^2,1)];
Y = C*X;
n = 2;
s = 2;
arr2avi('vid_example.mp4',Y,l,l)

% apply sid_dmd
[Ahat,Chat,Psi,Lam] = sid_dmd(Y,n,s);

% plot the results
imwrite(arr2img(real(Psi(:,1)),l,l)/0.02,'img_smode_1.png');
imwrite(arr2img(imag(Psi(:,1)),l,l)/0.02,'img_smode_2.png');
fg=figure;
hold on; box on; grid on;
plot(real(Lam(1).^(0:N)));
plot(imag(Lam(1).^(0:N)),'--');
xlabel('$k$','interpreter','latex')
ylabel('$\Re[\psi^k]$, $\Im[\psi^k]$','interpreter','latex')
xlim([0,N]);
saveas(fg,'img_tmode.eps');


function arr2avi(fpath,arr,nx,ny)
    v = VideoWriter(fpath,'MPEG-4');
    open(v);
    for j = 1:size(arr,2)
        writeVideo(v,arr2img(arr(:,j),nx,ny));
    end
    close(v);
end

function img = arr2img(arr,nx,ny)
    arr = reshape(arr,ny,nx);
    img = zeros(ny,nx,3);
    img(:,:,1) =-max(min(arr,0),-1)+min(max(arr,0), 1);
    img(:,:,2) = min(max(arr,0), 1);
    img(:,:,3) = min(max(arr,0), 1);
end