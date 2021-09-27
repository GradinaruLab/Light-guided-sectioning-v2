function [AmovingImage,AfixedImage]=adjust_brightness(movingImage,fixedImage,brightness1,brightness2,show_brightness)
% Anat Kahan, Cell Reports , 2021
%% brightness adjustment  for ' Register_2D_Zstack_2P_v2 '


 [AmovingImage]=adjust_image_brightness(movingImage,brightness1);
%check
 [AfixedImage]=adjust_image_brightness(fixedImage,brightness2);

 if show_brightness
figure
subplot (2,1,1), montage({AmovingImage,movingImage}), title('brightness original vs moving Image - GRIN')

subplot (2,1,2), montage({AfixedImage,fixedImage}), title('brightness original vs moving Image - GRIN')

 end