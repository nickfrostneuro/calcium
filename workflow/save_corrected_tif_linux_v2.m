function save_corrected_tif_linux_v2(correct_img, file, kk, dir_data, fr)
%save_correct_tif subfunction

     correct_img=uint16(correct_img);
     temp_dir =  dir_data{kk}; %convert cell to char
     dd = strcat(dir_data, '/xy_corrected/');
     ddchar = dd{kk};          %convert cell to char

       
     e = strcat(ddchar, 'frame_', num2str(fr,'%06d'), '_corrected', '.tif');
     imwrite (correct_img(:,:), e);
    