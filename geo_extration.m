%% % GEOMETRY OF THE EXTRACTION REGION
clear;
clc;

%select the plate
platenum = 3;
v = [0 1 15000];
x_ini = [0.0 1 0.0]; % Initial x coordinate of each plate
wire_size = 0.003; %cm 20% larger wires...
aper = 0.05; %cm 20% larger aperture
reg_cell = [1 0.5 1];
num_cell = round(reg_cell(platenum)/(wire_size + aper));
y_ini = [-1.5 1.5 3.5];
% Print the coordinates in Superfish format.
[fileID, errmsg] = fopen('input005.am','wt');
    formatSpec_header = ['&reg mat=0,voltage=' num2str(v(platenum)) ',ibound=-1 & \n'];
    formatSpec = '&po x=%6.4f,y=%6.4f & \n';
    
% File header and namelist for the superfish input
    fprintf(fileID,'DRATF Electrostatic Problem, Extraction Electrodes \n');
    fprintf(fileID,'&reg kprob=0,    ! Poisson or Pandira problem \n');
    fprintf(fileID,'xjfact=0.0,      ! Electrostatic problem \n');
    fprintf(fileID,'dx=0.004&       ! Mesh interval \n\n');
    fprintf(fileID, '&po x=0.0,y=0.0 &       ; Start of the boundary points \n');
    fprintf(fileID, '&po nt=5,radius=7.5,x=7.5,y=7.5 & \n');
    fprintf(fileID, '&po nt=5,radius=7.5,x=15,y=0.0 & \n');
    fprintf(fileID, '&po x=0.0,y=0.0 & \n\n');
for k=1:1;    
    fprintf(fileID, ['; ---------- PLATE ' num2str(k) '\n']);
    [CELL] = create_mesh_points(x_ini(k),num_cell, wire_size, aper);
for i=1:num_cell
    fprintf(fileID,formatSpec_header);
    %fwrite(fileID,formatSpec_header);
    for j=1:5
        % Swap x and y for GPT format, 1st and 2nd column in Coords.cell
        fprintf(fileID, formatSpec,CELL(i).Coords_cell(j,2)+x_ini(platenum),CELL(i).Coords_cell(j,1)+y_ini(platenum));
        %fwrite(fileID, ['&po x=' num2str(CELL(i).Coords_cell(j,1)) ' , y = ' num2str(CELL(i).Coords_cell(j,2)) ' & \n'], 'double');
    end
end

end

fclose(fileID);
