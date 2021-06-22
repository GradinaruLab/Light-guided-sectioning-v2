
function [theta_recovered, scale_recovered,tx,ty]= LGS_get_angle_and_scale(mytform)
% finds the angle and the scale of 2D registration using 'fitgeotrans' in
% 'Register_2D_Zstack_2P' 
% based on https://www.mathworks.com/matlabcentral/answers/167238-how-to-use-functions-fitgeotrans-and-transformpointsforward-instead-of-cp2tform-and-tformfwd

 t=mytform;
u = [0 1]; 
v = [0 0]; 
projectiveObj = projective2d(t.T);
[x, y] = transformPointsForward(projectiveObj, u, v);
dx = x(2) - x(1); 
dy = y(2) - y(1); 
angle = (180/pi) * atan2(dy, dx);
scale = 1 / sqrt(dx^2 + dy^2);

%% solution based on 'https://www.mathworks.com/help/images/find-image-rotation-and-scale.html?searchHighlight=fitgeotrans&s_tid=doc_srchtitle'
% that might fit better to fitgeotrans
tformInv = invert(t);
Tinv = tformInv.T;
ss = Tinv(2,1);
sc = Tinv(1,1);
ty=Tinv(3,2);
tx=Tinv(3,1);
scale_recovered = sqrt(ss*ss + sc*sc);
theta_recovered = atan2(ss,sc)*180/pi;

if abs(angle-theta_recovered)>0.001 
    disp('check reovered tform angle')
end
if abs(scale-scale_recovered)>0.001
     disp('check reovered tform scale')
end

return
