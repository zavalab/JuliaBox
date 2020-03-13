% Written by Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (victor.zavala@wisc.edu)

close all
clear
clc

% Importing data
n=3; s=20;
cstudy('data/dmmp.avi','lc_output/lc_dmmp',3,20)
cstudy('data/water.avi','lc_output/lc_water',3,20)

function cstudy(datapath,savepath,n,s)

    % obtain the data
    [Y,nx,ny] = avi2arr(datapath); Y=Y(:,1:71);
    fps = 30; trange = 0:size(Y,2)-1; dt = 1/fps;

    % identify the model
    [A,C,Psi,Lam] = sid_dmd(Y,n,s);

    % rearrange the spatiotemporal modes
    [~,p]=sort(abs(imag(Lam))-1e-6*real(Lam));
    Psi = Psi(:,p);
    Lam = Lam(p);

    % plot the results
    for i=1:n
        if imag(Lam(i))<0
            continue
        end
        imwrite(arr2img(real(Psi(:,i)),nx,ny),sprintf('%s_smode_%i.png',savepath,i))
        close all; fg=figure(1); fg.Position(3:4)= [300 100]; hold on; grid on; box on;
        fg.PaperPositionMode = 'auto';fig_pos = fg.PaperPosition;
        fg.PaperSize = [fig_pos(3) fig_pos(4)];
        plot(trange*dt,real(Lam(i).^trange),'k-'); xlim([trange(1) trange(end)]*dt);
        if Lam(i)~=0
            plot(trange*dt,imag(Lam(i).^trange),'k--'); xlim([trange(1) trange(end)]*dt);
        end
        saveas(fg,sprintf('%s_tmode_%i.pdf',savepath,i),'pdf');
        if imag(Lam(i))>0
            imwrite(arr2img(imag(Psi(:,i)),nx,ny),sprintf('%s_smode_%i.png',savepath,i+1))
            close all; fg=figure(1); fg.Position(3:4)= [300 100]; hold on; grid on; box on;
            fg.PaperPositionMode = 'auto';fig_pos = fg.PaperPosition;
            fg.PaperSize = [fig_pos(3) fig_pos(4)];
            plot(trange*dt,real(Lam(i).^trange),'k-');
            plot(trange*dt,imag(Lam(i).^trange),'k--'); xlim([trange(1) trange(end)]*dt);
            xlim([trange(1) trange(end)]*dt);
            saveas(fg,sprintf('%s_tmode_%i.pdf',savepath,i+1),'pdf');
        end
    end
end

function [X,Width,Height]=avi2arr(filepath)
    vd = VideoReader(filepath);
    Width = vd.Width; Height=vd.Height;
    X=[];
    while hasFrame(vd)
        X = [X reshape(double(rgb2gray(readFrame(vd)))/256,[],1)];
    end
end

function img = arr2img(arr,nx,ny)
    mag = norm(arr,Inf);
    arr = reshape(arr,ny,nx)/mag*sign(mean(arr(:)));
    img = zeros(ny,nx,3);
    img(:,:,1) =-max(min(arr,0),-1) + min(max(arr,0), 1);
    img(:,:,2) = min(max(arr,0), 1);
    img(:,:,3) = min(max(arr,0), 1);
end