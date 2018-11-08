fnames = {'hw2_3_1_a_50N40P10Z','hw2_3_1_a_80N30P30Z','hw2_3_1_c_both', ...
    'hw2_3_1_c_kzp','hw2_3_1_c_vmax','hw2_3_2_a_15','hw2_3_2_a_45','hw2_3_2_a_65','hw2_3_2_a_75'};
for i = 1:length(fnames)
    nnames = strcat(fnames{i},'.fig');
    open(nnames)
    print('-djpeg', fnames{i})
    close all
end