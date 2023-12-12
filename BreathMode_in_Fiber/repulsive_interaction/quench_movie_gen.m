load wavefunctions_t.mat;
[M,~,~] = size(Pw);
for ind = 1:M
    imagesc(squeeze(Pw(ind,:,:)));axis equal
    %zlim([0 0.08])
    shading interp
    F(ind) = getframe(gcf);
end
%drawnow
video_title = strcat('triangle_movie');
writerObj = VideoWriter(video_title);
writerObj.FrameRate = 10;
open(writerObj);
writeVideo(writerObj,F(2:length(F)));
close(writerObj);