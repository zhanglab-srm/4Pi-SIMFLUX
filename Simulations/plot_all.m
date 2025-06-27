close all
clear
clc

data = load('Fitting_std_single_len.mat');
figure
z = data.data_save.z_save;
x_pres = data.data_save.STD_x;
y_pres = data.data_save.STD_y;
z_pres = data.data_save.STD_z;
plot(z, x_pres, z, y_pres, z, z_pres);
legend('x', 'y', 'z')
title('3D-SMLM localization precision')
axis([-500 500 0 20])

data = load('Fitting_std_4pi.mat');
figure
z = data.data_save.z_save;
x_pres = data.data_save.STD_x;
y_pres = data.data_save.STD_y;
z_pres = data.data_save.STD_z;
plot(z, x_pres, z, y_pres, z, z_pres);
legend('x', 'y', 'z')
title('4Pi-SMLM localization precision')
axis([-500 500 0 20])

data = load('Fitting_std_4pi_SIMFLUX.mat');
figure
z = data.data_save.z_save;
x_pres = data.data_save.STD_x;
y_pres = data.data_save.STD_y;
z_pres = data.data_save.STD_z;
plot(z, x_pres, z, y_pres, z, z_pres);
legend('x', 'y', 'z')
title('4Pi-SIMFLUX localization precision')
axis([-500 500 0 20])

data = load('Fitting_std_4pi_SIMFLUX_compare.mat');
figure
z = data.data_save.photons;
x_pres = data.data_save.STD_x;
y_pres = data.data_save.STD_y;
z_pres = data.data_save.STD_z;
plot(z, x_pres)
loglog(z, x_pres)
hold on
plot(z, y_pres)
loglog(z, y_pres)
hold on
plot(z, z_pres)
loglog(z, z_pres)
legend('x', 'y', 'z')
title('4Pi-SIMFLUX localization precision vs photons')