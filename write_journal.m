filename=['result_' datestr(today,'yyyymmdd') '.txt'];

fileID=fopen(filename,'a');

flist = dir('pre_setup.m')

copyfile(flist(1).name,filename);

fprintf(fileID,'\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nResults:\n\n');

fprintf(fileID,'Outlet temperature\n %.2f\n\n',t_out);

fprintf(fileID,'Thermal efficiency\n %.4f\n\n',thermal_eff);

fprintf(fileID,'Electric efficiency\n %.4f\n\n',ele_eff);

fprintf(fileID,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');