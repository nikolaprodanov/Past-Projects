clear all;
close all;
clc
!--------------- File Making -------------------
element='Si'; !Modify if the simulation is made with a different element
!-----------------------------------------------
prefix='OUTPUT_DOS_of_';
prefix_2='OUTPUT_electron_spectrum_';
suffix='.dat';
!------------- Adding the files ----------------
file=strcat(prefix,element,suffix);
file_2=strcat(prefix_2,element,suffix);
!--------------- Reading the Files -------------
A=dlmread('OUTPUT_total_all.dat','',2,0); !the number two means that it is skipping the first two lines!
B=dlmread( file ,'',1,0);
C=dlmread( file_2 ,'',2,0);
!--------------- Assigning variables -----------
energy=B(:,1);
DOS=B(:,2);
k_vector=B(:,3);
m_eff=B(:,4);

time=A(:,1);
energy_total=A(:,12);
energy_photons=A(:,6)
energy_electrons=A(:,7);
energy_holes_kin=A(:,8);
energy_holes_pot=A(:,9);
energy_positrons=A(:,10);
energy_atoms=A(:,11);

num_photons=A(:,2);
num_electrons=A(:,3);
num_holes=A(:,4);
num_positrons=A(:,5);

time_2=C(:,1);
energy_2=C(:,2);
distribution=C(:,3);
!--------------- Plotting Figures ---------------
!--------------- DOS figures --------------------
figure
plot(energy,DOS,'linewidth',2)
title('Density of Electron States')
set(gca,'fontsize',14,'linewidth',1)
xlabel('Energy [eV]');
ylabel('DOS [eV^-1]');
saveas(gcf,'OUTPUT_DOS','png');
close

figure
plot(energy,k_vector,'linewidth',2)
title('k vector vs Energy')
set(gca,'fontsize',14,'linewidth',1)
xlabel('Energy [eV]');
ylabel('k-vector [m^-1]');
saveas(gcf,'OUTPUT_DOS_k_vector','png');
close

figure
plot(energy,m_eff,'linewidth',2)
title('Effective mass vs Energy')
set(gca,'fontsize',14,'linewidth',1)
xlabel('Energy [eV]');
ylabel('Effective mass [me]');
ylim([0 10])
saveas(gcf,'OUTPUT_DOS_effective_mass','png');
close 

!------------------- Total Energy ------------------
figure
plot(time,energy_total,'linewidth',2)
hold on
plot(time,energy_photons,'linewidth',2)
plot(time,energy_electrons,'linewidth',2)
plot(time,energy_holes_kin,'linewidth',2)
plot(time,energy_holes_pot,'linewidth',2)
plot(time,energy_positrons,'linewidth',2)
plot(time,energy_atoms,'linewidth',2)
hold off
legend({'Total','Photons','Electrons','Holes (kin)','Holes (pot)','Positrons','Atoms'})
title('Time evolution of particles energy')
set(gca,'fontsize',14,'linewidth',1)
xlabel('time [fs]');
ylabel('Energy [eV]');
saveas(gcf,'OUTPUT_total_energies','png');
close

!------------------- Number of particles -----------
figure
plot(time,num_photons,'linewidth',2)
hold on
plot(time,num_electrons,'linewidth',2)
plot(time,num_holes+0.1,'linewidth',2)
plot(time,num_positrons,'linewidth',2)
hold off
legend({'Photons','Electrons','Holes','Positrons'})
title('Time evolution of particle count')
set(gca,'fontsize',14,'linewidth',1)
xlabel('time [fs]');
ylabel('Number of particles');
saveas(gcf,'OUTPUT_total_numbers','png');
close

!------------------- Energy of Electrons Distribution --------
filename = 'OUTPUT_electron_spectrum_1d_Z.dat';

data = readtable(filename, 'HeaderLines', 3);

dataArray = table2array(data);

time_to_read = 100;

filteredData = dataArray(dataArray(:, 1) == time_to_read, :);

figure;
plot(filteredData(:, 2), filteredData(:, 3), 'LineWidth', 2);
title(['External electron Energy Distribution at ',num2str(time_to_read), ' fs'])
set(gca,'fontsize',14,'linewidth',1)
xticks(0:10:90)
!ylim([0 0.0006])
xlabel('Energy [eV]');
ylabel('Distribution [eV^-1]');
saveas(gcf,'OUTPUT_electron_spectrum','png');
close
