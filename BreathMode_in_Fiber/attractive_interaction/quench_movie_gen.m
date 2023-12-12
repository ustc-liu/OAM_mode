load wavefunctions_t.mat;
[M,~,~] = size(Pw);
for ind = 1:M
    surf(squeeze(Pw(ind,:,:)));%axis equal
    zlim([0 0.06])
    shading interp
    F(ind) = getframe(gcf);
end
%drawnow
video_title = strcat('disk_movie_2');
writerObj = VideoWriter(video_title);
writerObj.FrameRate = 10;
open(writerObj);
writeVideo(writerObj,F(2:length(F)));
close(writerObj);