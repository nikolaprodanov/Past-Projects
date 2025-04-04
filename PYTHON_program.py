# ************************************************************************
#                            INPUT DATA 
# ************************************************************************

# The variable we want to vary for the study of the electron yield: (X-AXIS)
chosen_variable = 1 
# 0 = material thickness(z2)
# 1 = incoming particle irradiation energy 
# 2 = polar angle of incidence
# 3 = Monte Carlo iterations

# Ranges for the study of the Electron Spectrum in units of the chosen variable
variable_min, variable_max = 0, 1e4

# Time to read the electron spectrum [fs]
time_to_read = 15

# Material Properties
work_function = 4.67 # [eV] work function of the material. Needed for the energy shift in spectrum.

# ************************************************************************
#                       E-YIELD Analysis Settings 
# ************************************************************************

# Noise filtering settings: Savitsky Golay filter
noise_filtering = True
polinomial_order_filter = 2
window_length_filter = 20
plot_with_noise_error = True

# Time settings
err_time = 0.001 # Fixed Variable: Error in time reading. Must be smaller than the dt_step in TREKIS-4 which is usually 1.0 fs

# End time of cascades
f = 0.95 # [1] for the Threshold Method
treshold_angle = 10 # [degrees] for the Derivative Threshold Method
normalization = False


# ************************************************************************
#                       E-YIELD production settings 
# ************************************************************************

# Choose which plots & Analysis to make
make_plot_e_spectrum = True
make_plot_e_yield = True
make_data_and_plot_end_of_cascades = True
make_data_and_plot_end_of_eyield = True
make_GIFs = False

# Images output settings
keep_e_spectrum_image = False # not used
keep_e_yield_image = False # not used 
keep_cascade_time_vs_variable = True
keep_electrons_vs_time = False
keep_eyield_saturation_vs_variable = True
keep_eyield_vs_time = False
keep_GIF_generated_images = False

# Print Folders settings
print_folders_and_significant_data = True

# Styling settings (purely visual)
Styling = True

# Which version of TREKIS-4 are we using
fedor_version = False # Set it to True only if you are comparing older versions of TREKIS-4

# ************************************************************************
#                           Python LibrariesS 
# ************************************************************************

import os 
import sys
import scipy
from scipy import signal
from scipy.optimize import curve_fit 
import numpy as np 
import matplotlib.pyplot as plt 
import time 
import imageio
import math

# ************************************************************************
#                           FILE NAMES
# ************************************************************************

file_e_spectrum = 'OUTPUT_electron_spectrum_1d_Z.dat'
file_e_density = 'OUTPUT_electron_density_1d_Z_Si.dat' # edit this so it becomes for any material possible to work with
file_OUTPUT_total = 'OUTPUT_total_all.dat'
directory_data_folders = 'TREKIS_OUTPUTS'

# ************************************************************************
#                       Definitions of Functions 
# ************************************************************************

# Extracting Data 
def extract_data_from_spectrums():
    data = [] # variable, e_yield, e_energy, e_spectrum, time, d, theta, barrier, material, irradiation energy, error e-yield
    for subfolder in os.listdir(directory_data_folders):
        subfolder_path = os.path.join(directory_data_folders, subfolder)
        
        if os.path.isdir(subfolder_path):

            # Reading the current variable and other settings
            current_variable = read_variable(chosen_variable,subfolder_path) 
            virt_d = read_INPUT_DATA_material_z2(subfolder_path)
            virt_theta = read_INPUT_DATA_theta(subfolder_path)
            virt_barrier = read_INPUT_DATA_barrier(subfolder_path) 
            virt_material = read_INPUT_DATA_material(subfolder_path)
            virt_irr_energy = read_INPUT_DATA_irradiation_energy(subfolder_path)

            # ***************************************************************
            # Reading the electron spectrum and energy

            with open(os.path.join(subfolder_path, file_e_spectrum), 'r') as file:
                lines = file.readlines()
            
            current_time = None
            e_energy = []
            e_spectrum = []

            for line in lines:
                if line.startswith('#') or line.strip() == '':
                    continue # This skips commented or empty lines
                
                parts = line.split()
                if parts[2] != 'Infinity': # if an error occurred in TREKIS
                    s_time, s_energy, s_spectrum = float(parts[0]), (float(parts[1]) - work_function), float(parts[2]) # s stands for single, CHANGE SPECTRUM COLUMN HERE
                else:
                    s_time, s_energy, s_spectrum = float(parts[0]), float(0), float(parts[2])
                if s_energy < 0:
                    continue # Skips the shifted energies in the electron spectrum

                if current_time is None:
                    current_time = s_time

                if s_time != current_time:
                    
                    # noise filter to spectrum # calculate the yield
                    if noise_filtering == True:
                        e_yield_noisy = scipy.integrate.simpson(e_spectrum, e_energy)
                        e_spectrum = signal.savgol_filter(e_spectrum, window_length=window_length_filter, polyorder=polinomial_order_filter, mode="nearest")
                    
                    # Calculating the Electron Yield with the Simpson Integral
                    e_yield = scipy.integrate.simpson(e_spectrum, e_energy)
                    if noise_filtering == True:
                        error_e_yield = abs(e_yield_noisy-e_yield)
                    else:
                        error_e_yield = 0

                    # Calculating the noisy deviation in the electron yield
                    

                    # save the data
                    # variable, e_yield, e_energy, e_spectrum, time, d, theta, barrier, material, irradiation energy, error e-yield
                    data.append((current_variable, e_yield, e_energy, e_spectrum, current_time, virt_d, virt_theta, virt_barrier, virt_material, virt_irr_energy, error_e_yield ))

                    # reset the values to none for other times
                    current_time = s_time
                    e_energy = []
                    e_spectrum = []

                e_energy.append(s_energy)
                e_spectrum.append(s_spectrum)
            
            # Append the very last set of data, the very last line for the last time:
            if s_energy:
                data.append((current_variable, e_yield, e_energy, e_spectrum, current_time, virt_d, virt_theta, virt_barrier, virt_material, virt_irr_energy, error_e_yield))
    #radiation_energy = read_INPUT_DATA_irradiation_energy(subfolder_path)
    return data

def extract_data_from_electron_density():
    data = [] # variable, time, Z-coordinate, density
    for subfolder in os.listdir(directory_data_folders):
        subfolder_path = os.path.join(directory_data_folders, subfolder)
        
        if os.path.isdir(subfolder_path):

            # Reading the current variable and other settings
            current_variable = read_variable(chosen_variable,subfolder_path) 

            # ***************************************************************
            # Reading the electron density in z coordinates

            with open(os.path.join(subfolder_path, file_e_density), 'r') as file:
                lines = file.readlines()
            
            current_time = None
            coordinates = []
            densities = []

            for line in lines:
                if line.startswith('#') or line.strip() == '':
                    continue # This skips commented or empty lines
                
                parts = line.split()
                s_time, s_coordinate, s_density = float(parts[0]), (float(parts[1])), float(parts[2]) # s stands for single, CHANGE SPECTRUM COLUMN HERE

                if s_coordinate < 0:
                    continue # Skips the negative values of Z & zero, due to stamping error, fix this in future

                if current_time is None:
                    current_time = s_time

                if s_time != current_time:
                    
                    # noise filter to spectrum # calculate the yield
                    #if noise_filtering == True:
                    #    e_yield_noisy = scipy.integrate.simpson(e_spectrum, e_energy)
                    #    e_spectrum = signal.savgol_filter(e_spectrum, window_length=window_length_filter, polyorder=polinomial_order_filter, mode="nearest")
                    
                    # Calculating the Electron Yield with the Simpson Integral
                    #e_yield = scipy.integrate.simpson(e_spectrum, e_energy)
                    #if noise_filtering == True:
                    #    error_e_yield = abs(e_yield_noisy-e_yield)
                    #else:
                    #    error_e_yield = 0

                    # Calculating the noisy deviation in the electron yield
                    

                    # save the data
                    # variable, e_yield, e_energy, e_spectrum, time, d, theta, barrier, material, irradiation energy, error e-yield
                    data.append((current_variable, current_time, coordinates, densities ))

                    # reset the values to none for other times
                    current_time = s_time
                    coordinates = []
                    densities = []

                coordinates.append(s_coordinate)
                densities.append(s_density)
            
            # Append the very last set of data, the very last line for the last time:
            if s_coordinate:
                data.append((current_variable, current_time, coordinates, densities ))
    #radiation_energy = read_INPUT_DATA_irradiation_energy(subfolder_path)
    return data

# Read Variables
def read_variable(integer, file_path):
    # 0 corresponds to the material thickness (that is Z2 border of the material)
    # 1 corresponds to the irradiation energy
    # 2 correpsonds to the angle theta of incidence of the incoming particle
    # 3 corresponds to the number of monte carlo iterations
    if integer == 0:
        return float(read_INPUT_DATA_material_z2(file_path))
    if integer == 1:
        return float(read_INPUT_DATA_irradiation_energy(file_path))
    if integer == 2:
        return float(read_INPUT_DATA_theta(file_path))
    if integer == 3:
        return float(read_NUMERICAL_PARAMETERS_MC_iterations(file_path))
    else:
        if integer != 0 and integer !=1 and integer !=2:
            print(f'{prfx}Error: The integer you inserted is wrong! you can insert 0,1,2.')
            exit()
        else:
            print(f'{prfx}Error: file_path or something else')
            exit()

def variable_settings(integer):
    # 0 corresponds to the material thickness (that is Z2 border of the material)
    # 1 corresponds to the irradiation energy
    # 2 correpsonds to the angle theta of incidence of the incoming particle
    # 3 corresponds to the number of monte carlo iterations
    # returns 'name', 'label', 'unit'
    if integer == 0:
        return 'film thickness', 'd', 'A'
    if integer == 1:
        return 'irradiation energy', 'E_prtcl', 'eV'
    if integer == 2:
        return 'polar angle of incidence', 'theta', ''
    if integer == 3:
        return 'Monte Carlo Iterations', 'MC_iter', ''
    else:
        if integer != 0 and integer !=1 and integer !=2 and integer !=3:
            print(f'{prfx}Error: The integer you inserted is wrong! you can insert 0,1,2,3.')
            exit()
        else:
            print(f'{prfx}Error: unknown')
            exit()

# Read files
def read_INPUT_DATA_irradiation_energy(file_path):
    #Reads from INPUT_DATA the irradiation energy of the incoming particle
    # number of line counting from 0 in TREKIS-4 is 21: fedors version is 22
    num_line = 21
    if fedor_version == True:
        num_line = 22
    try:
        with open(file_path+'\INPUT_DATA.txt', 'r') as f:
            lines = f.readlines()
            if len(lines) < 22:
                myprint(f'{prfx}Warning: {file_path} is constructed wrong.')
                return None
            irradiation_energy_line = lines[num_line].strip()
            irradiation_energy = irradiation_energy_line.split()[0]
            return float(irradiation_energy)
    except ValueError:
        print(f'{prfx}Error: Line {num_line} in {file_path} does not contain a valid number.')
    except Exception as e:
        print(f'{prfx}Error reading {file_path}: {e}')
    return None

def read_INPUT_DATA_material_z2(file_path):
    #Reads from INPUT_DATA the second variable z2 that indicates the end of the material
    # number of line counting from 0 in TREKIS-4 is 9: fedors version is also
    num_line = 9
    try:
        with open(file_path+'\INPUT_DATA.txt', 'r') as f:
            lines = f.readlines()
            if len(lines) < 22:
                myprint(f'{prfx}Warning: {file_path} is constructed wrong.')
                return None
            z2_line = lines[num_line].strip()
            z2 = z2_line.split()[1]
            return float(z2)
    except ValueError:
        print(f'{prfx}Error: Line {num_line} in {file_path} does not contain a valid number.')
    except Exception as e:
        print(f'{prfx}Error reading {file_path}: {e}')
    return None

def read_INPUT_DATA_material(file_path):
    #Reads from INPUT_DATA the material
    # number of line counting from 0 in TREKIS-4 is 1: also in fedors version
    num_line = 1
    try:
        with open(file_path+'\INPUT_DATA.txt', 'r') as f:
            lines = f.readlines()
            if len(lines) < 22:
                myprint(f'{prfx}Warning: {file_path} is constructed wrong.')
                return None
            material_line = lines[num_line].strip()
            material = material_line.split()[0]
            return (material)
    except ValueError:
        print(f'{prfx}Error: Line {num_line} in {file_path} does not contain a valid number.')
    except Exception as e:
        print(f'{prfx}Error reading {file_path}: {e}')
    return None

def read_INPUT_DATA_theta(file_path):

    #Reads from INPUT_DATA the angle of incidence of the incoming particle
    # number of line counting from 0 in TREKIS-4 is 19: in fedors version it is 20
    num_line = 19
    if fedor_version == True:
        num_line = 20
    try:
        with open(file_path+'\INPUT_DATA.txt', 'r') as f:
            lines = f.readlines()
            if len(lines) < 22:
                myprint(f'{prfx}Warning: {file_path} is constructed wrong.')
                return None
            material_line = lines[num_line].strip()
            material = material_line.split()[0]
            return (material)
    except ValueError:
        print(f'{prfx}Error: Line {num_line} in {file_path} does not contain a valid number.')
    except Exception as e:
        print(f'{prfx}Error reading {file_path}: {e}')
    return None

def read_INPUT_DATA_barrier(file_path):
    #Reads from INPUT_DATA the potential barrier of the solid
    # number of line counting from 0 in TREKIS-4 is 13: in fedors version it is non-existent so we will set it to step
    num_line = 13
    if fedor_version == True:
        return 'step'
    try:
        with open(file_path+'\INPUT_DATA.txt', 'r') as f:
            lines = f.readlines()
            if len(lines) < 22:
                myprint(f'{prfx}Warning: {file_path} is constructed wrong.')
                return None
            barrier_line = lines[num_line].strip()
            barrier = int(barrier_line.split()[0])
            if barrier == 0:
                return 'step'
            if barrier == 1:
                return 'Eckart'
            else:
                print(f'{prfx}The barrier red is {barrier}. It should be 0 or 1.')
                exit()
    except ValueError:
        print(f'{prfx}Error: Line {num_line} in {file_path} does not contain a valid number.')
    except Exception as e:
        print(f'{prfx}Error reading {file_path}: {e}')
    return None

def read_NUMERICAL_PARAMETERS_MC_iterations(file_path):
    #Reads from NUMERICAL_PARAMETERS the number of montecarlo iterations
    # number of line counting from 0 in TREKIS-4 is 13: in fedors version it is non-existent so we will set it to step
    num_line = 1
    try:
        with open(rf'{file_path}\NUMERICAL_PARAMETERS.txt', 'r') as f:
            lines = f.readlines()
            if len(lines) < 22:
                myprint(f'{prfx}Warning: {file_path} is constructed wrong.')
                return None
            separate_line = lines[num_line].strip()
            MC_iterations = int(separate_line.split()[0])
            return MC_iterations
    except ValueError:
        print(f'{prfx}Error: Line {num_line} in {file_path} does not contain a valid number.')
    except Exception as e:
        print(f'{prfx}Error reading {file_path}: {e}')
    return None

def get_max_time(path, filename):
    file_path = os.path.join(path, filename)
    max_time = float('-inf')
    try:
        with open(file_path, 'r') as file:
            #Skips the first three lines since they are commented
            for _ in range(3):
                next(file)
            #Strips the lines
            for line in file:
                if not line.strip():
                    continue
                # splits the lines 
                columns = line.split()
                if len(columns) < 2:
                    continue
                #Get the time from the file
                time = float(columns[0])
                if max_time == float('-inf'):
                    max_time = time
                if time > max_time:
                    max_time = time
        return max_time
    except Exception as e:
        print(f'Error reading file {file_path}: {e}')

def read_OUTPUT_total_all(path, filename, costum_column):
    file_path = os.path.join(path, filename)
    x = []
    y = []
    try:
        with open(file_path, 'r') as file:
            #Skips the first three lines since they are commented
            for _ in range(2):
                next(file)
            #Strips the lines
            for line in file:
                if not line.strip():
                    continue
                # splits the lines 
                columns = line.split()
                if len(columns) < 2:
                    continue
                x.append(float(columns[0]))
                y.append(float(columns[costum_column]))
        return x,y
    except Exception as e:
        print(f'Error reading file {file_path}: {e}')

# Checking if the folders and files are compatible
def check_if_directory_exists(directory_data_folders):
    if not os.path.isdir(directory_data_folders):
        print(f'{prfx}The directory {directory_data_folders} does not exist.')
        print(f'{prfx}Exiting')
        sys.exit(1)

def count_the_subfolders(directory_data_folders):
    i = 0
    for subfolder in os.listdir(directory_data_folders):
        subfolder_path = os.path.join(directory_data_folders, subfolder)
        if os.path.isdir(subfolder_path):
            i = i + 1
    return i

# Saving Data, Plotting and GIFs

def plot_electron_spectrum(data, target_time):

    if make_plot_e_spectrum == False:
        return

    print(f'{prfx}Making Plot: Electron Spectrum with varying {variable_name} ({variable_unit}) = [{variable_min},{variable_max}]')

    for i in range(0, len(data) - 1):
        var = data[i][0]
        time = data[i][4]
        if var >= variable_min and var <= variable_max and abs(time - target_time) <= err_time:
            
            e_energy = data[i][2]
            e_spectrum = data[i][3]
            material = data[i][8]

            # -------------------------------------------------------------- Sub function: Normalize to the maximum in some region
            #g_energy = []
            #g_spectrum = []
            #for j in range(0,len(data[i][3]) - 1):
            #    if e_energy[j] <= 30:
            #        g_energy.append(e_energy[j])
            #        g_spectrum.append(e_spectrum[j])
            #    else:
            #        break
            #g_max_spectrum = max(g_spectrum)
            #for j in range(0, len(g_energy) - 1):
                #g_spectrum[j] = g_spectrum[j]/g_max_spectrum

            #plt.plot(g_energy, g_spectrum, color = 'r',label = 'TREKIS')
            # ------------------------------------------------------------- End of Sub function
            plt.plot(e_energy, e_spectrum, label = f'{variable_label} = {var} {variable_unit}')
            # Saving each spectrum in a different labeled file
            with open(f'e_spectrum_at_{time_to_read}_fs_{variable_label}={var}.dat', 'w') as file:
                for i in range(0, len(e_energy)):
                    file.write(f'{e_energy[i]} {e_spectrum[i]}\n')
    
    plt.legend()
    plt.xlim(0,100)
    plt.ylim(0)
    material = data[0][8]
    plt.title(f'Outer Electron Spectrum of {material} at {(time_to_read)}fs')
    plt.xlabel('Electron Energy (eV)')
    plt.ylabel('Number of electrons(1/inc.photon/eV)')
    plt.savefig(f'IMAGE_Outer_e_spectrum_of_{material}_for_{variable_label}_{variable_min}_{variable_max}.png')
    plt.close()

def plot_electron_yield(data, target_time):

    if make_plot_e_yield == False:
        return

    # data = [] variable, e_yield, e_energy, e_spectrum, time, d, theta, barrier, material, N_electrons, irradiation energy, error e-yield
    print(f'{prfx}Making Plot: Electron Yield VS {variable_name} ({variable_unit})')

    var = []
    e_yield = []
    error = []

    for i in range(0, len(data) - 1):
        time = data[i][4]
        if abs(time - target_time) <= err_time:

            var.append(data[i][0])
            e_yield.append(data[i][1])
            error.append(data[i][10])
    material = data[0][8]
    if plot_with_noise_error == True:
        plt.errorbar(var, e_yield, error, marker = 'o', linestyle = '-')
    else:    
        plt.plot(var , e_yield, marker = 'o', linestyle = '-')
    plt.title(f'Electron yield of {material} for different {variable_name} at {(time_to_read)} fs')
    plt.xlabel(f'{variable_name} ({variable_unit})')
    plt.ylabel('Electron yield(1/inc.photon)')
    plt.ylim(0)
    plt.savefig(f'IMAGE_e_yield_at_{time_to_read}_fs.png')
    plt.close()
    # Saving the data in a file
    with open(f'electron_yield_at_{time_to_read}_fs.dat', 'w') as file:
        for var,e_yield in zip(var, e_yield):
            file.write(f'{var} {e_yield}\n')

def make_GIF_electron_spectrum_evolution(data):

    if make_GIFs == False:
        return

    print(f'{prfx}Making GIF: Electron Spectrum')

    output_directory = 'IMAGES_e_spectrum_evolution'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    modified_data = data

    modified_data.sort(key = lambda x: x[4]) # sorts the data by time

    filenames = []

    fig, ax = plt.subplots()
    for i in range(0, len(modified_data) - 1):
        
        time_i = data[i][4]
        e_energy = data[i][2]
        e_spectrum = data[i][3]
        var = data[i][0]
        material = data[i][8]
        radiation_energy = data[i][9]

        # Plot
        ax.plot(e_energy, e_spectrum, label = f'{variable_label} = {var} {variable_unit}')

        # Check if the following time is different than the current and therefore save the current time frame
        time_next = data[i+1][4]
        if time_next != time_i:
            ax.legend()
            ax.set_xlim(0,radiation_energy)
            ax.set_ylim(0)
            ax.set_title(f'Outer Electron Spectrum of {material} at {(time_i)}fs')
            ax.set_xlabel('Electron Energy (eV)')
            ax.set_ylabel('Number of electrons(1/inc.photon/eV)')

            filename = os.path.join(output_directory, f'e_spectrum_{time_i}.png')
            filenames.append(filename)
            plt.savefig(filename)
            plt.close(fig) # close figure so for the next time step we have a new plot
            fig, ax = plt.subplots() # open new figure

    plt.close(fig) # there is one extra blank figure which needs to be closed

    # Create a GIF from the saved images
    with imageio.get_writer('GIF_e_spectrum_evolution.gif', mode = 'I', duration = 500) as writer:
        for filename in filenames:
            image = imageio.v2.imread(filename)
            writer.append_data(image)

    # Delete the images
    if keep_GIF_generated_images == False:
        for filename in filenames:
            os.remove(filename)
    
def make_GIF_electron_yield_evolution(data):

    if make_GIFs == False:
        return

    print(f'{prfx}Making GIF: Electron Yield')

    output_directory = 'IMAGES_e_yield_evolution'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    modified_data = data

    modified_data.sort(key = lambda x: x[4]) # sorts the data by time

    filenames = []

    e_yield = []
    var = []
    fig, ax = plt.subplots()
    for i in range(0, len(modified_data) - 1):
        
        time_i = data[i][4]
        e_yield.append(data[i][1])
        var.append(data[i][0])
        material = data[i][8]

        # Check if the following time is different than the current and therefore save the current time frame
        time_next = data[i+1][4]
        if time_next != time_i:
            # Plot
            ax.plot(var, e_yield, marker = 'o', linestyle = '-')
            ax.set_ylim(0, 0.01)
            material = data[i][8]
            ax.set_title(f'Electron yield of {material} for different {variable_name} at {time_i} fs')
            ax.set_xlabel(f'{variable_name} ({variable_unit})')
            ax.set_ylabel('Electron yield(1/inc.photon)')

            filename = os.path.join(output_directory, f'e_yield_{time_i}.png')
            filenames.append(filename)
            plt.savefig(filename)
            plt.close(fig)
            fig, ax = plt.subplots()
            e_yield = []
            var = []
            #close figure so for the next time step we have a new plot

    plt.close(fig) # there is one extra blank figure which needs to be closed

    # Create a GIF from the saved images
    with imageio.get_writer('GIF_e_yield_evolution.gif', mode = 'I', duration = 500) as writer:
        for filename in filenames:
            image = imageio.v2.imread(filename)
            writer.append_data(image)

    # Delete the images
    if keep_GIF_generated_images == False:
        for filename in filenames:
            os.remove(filename)

def data_and_plot_end_of_cascades_vs_variable(): # lists = ((var1,t_end_1), ...)

    if make_data_and_plot_end_of_cascades == False:
        return None

    print(f'{prfx}Making Analysis & Plot: Time cascades of electrons in the Simulation Box')

    dN_dt = []
    data_cascades = []

    # The program goes back in the folder to read the file: OUTPUT_total_all.dat
    for subfolder in os.listdir(directory_data_folders):
        subfolder_path = os.path.join(directory_data_folders, subfolder)
        if os.path.isdir(subfolder_path):

            # Readings for later plots
            material = read_INPUT_DATA_material(subfolder_path)
            current_variable = read_variable(chosen_variable, subfolder_path)
        
            # Get the number of electrons (N_e) at each time t
            t, N_e = read_OUTPUT_total_all(subfolder_path, file_OUTPUT_total, 2)

            # Normalization
            if normalization == True:
                N_e, N_e_transform = my_renormalization(N_e)
                t, t_transform = my_renormalization(t)

            # Calculating the derivative
            dN_dt = my_derivative(N_e, t)

            # Find the time cascade end
            t_cascade_end = None
            for i in range(len(N_e)):
                #if t_cascade_end == None and abs(math.atan(dN_dt[i])*180/math.pi) < treshold_angle: # DIFFERENT METHOD ADD IN FUNCTION
                if t_cascade_end == None and N_e[i] >= f * N_e[-1]:
                    t_cascade_end = t[i]
                    break
            
            # Plotting
            # normalizing them to each other:
            max_derivative = max(dN_dt)
            for i in range(len(N_e)):
                N_e[i] = N_e[i]/N_e[-1]
                dN_dt[i] = dN_dt[i]/max_derivative
            plt.plot(t, N_e, linestyle='-', label='$N_{e}$ '+ f'for {variable_label} = {current_variable} {variable_unit}')
            plt.plot(t, dN_dt, linestyle='-', label = 'd$N_{e}$/dt')
            plt.axvline(x = t_cascade_end, linestyle='--', label = 'End of Cascade Line')

            if normalization == True:
                plt.xlim(0,1)
                plt.ylim(-0.05,1.05)
                plt.title(f'Normalized evolution of electrons in the Simulation Box for {variable_label} = {current_variable} {variable_unit}')
                plt.xlabel('time - Normalized')
                plt.ylabel('Number of electrons - Normalized')
            else:
                #plt.xlim(0)
                plt.ylim(0)
                plt.title(f'Evolution of electrons in the Simulation Box for {variable_label} = {current_variable} {variable_unit}')
                plt.ylabel('Number of electrons')
                plt.xlabel('time (fs)')
            
            plt.legend()
            if keep_electrons_vs_time == True:
                plt.savefig(f'IMAGE_e_evolution_{variable_name}={current_variable}.png')
            plt.close()

            # Anti - Normalization of the data
            if normalization == True:
                N_e = my_anti_renormalization(N_e, N_e_transform)
                t = my_anti_renormalization(t, t_transform)
                dN_dt = my_anti_renormalization(dN_dt, t_transform/N_e_transform)
                if t_cascade_end != None:
                    t_cascade_end = single_variable_renorm_to_antirenorm(t_cascade_end, t_transform)
        
            if t_cascade_end != None:
                data_cascades.append((current_variable, t_cascade_end))
                #print(f'{prfx}For {variable_name} = {current_variable}{variable_unit} : Cascade end is at {t_cascade_end} fs')
            else:
                print(f'{prfx}For {variable_name} = {current_variable}{variable_unit} : Cascade has not converged.')
    
    # Saving the data in a file
    data_cascades.sort(key=lambda x: x[0])
    with open(f'cascade_time_vs_{variable_name}.dat', 'w') as file:
        for i in range(len(data_cascades)):
            file.write(f'{data_cascades[i][0]}      {data_cascades[i][1]}\n')

    return data_cascades

def data_and_plot_end_of_eyield_vs_variable(data): # lists = ((var1, t_end_1), ...)
    
    if make_data_and_plot_end_of_eyield == False:
        return None

    print(f'{prfx}Making Analysis & Plot: Electron Yield Saturation')

    # Sort for Algorithm
    data.sort(key=lambda x: x[4]) # sort data by time
    data.sort(key=lambda x: x[0]) # sort by variable
    
    data_eyield_end = []
    eyield = []
    t = []
    deyield_dt = []

    for i in range(0, len(data) - 1):

        t_i = data[i][4]
        eyield_i = data[i][1]
        var_i = data[i][0]
        var_next = data[i+1][0]

        eyield.append(eyield_i)
        t.append(t_i)

        if var_next != var_i or i == len(data) - 2:
            # its time to analyze all the data for a given variable
            if normalization == True:
                eyield, eyield_transform = my_renormalization(eyield)
                t, t_transform = my_renormalization(t)

            # Calculating the derivative
            deyield_dt = my_derivative(eyield, t)

            # Find the time cascade end
            t_eyield_end = None
            for i in range(len(eyield)):
                if t_eyield_end == None and eyield[i] >= f * eyield[-1]:
                    t_eyield_end = t[i]
                    break
            
            # Plotting
            plt.plot(t, eyield, linestyle='-', label='Yield')
            plt.plot(t, deyield_dt, linestyle='-', label = 'dY/dt')
            plt.axvline(x = t_eyield_end, linestyle='--', label = 'Saturation')

            if normalization == True:
                plt.xlim(0,1)
                plt.ylim(-0.05,1.05)
                plt.title(f'Normalized evolution of the electron Yield for {variable_label} = {var_i} {variable_unit}')
                plt.xlabel('time - Normalized')
                plt.ylabel('Yield - Normalized')
            else:
                plt.xlim(0)
                plt.ylim(0)
                plt.title(f'Evolution of the electron Yield for {variable_label} = {var_i} {variable_unit}')
                plt.ylabel('Yield')
                plt.xlabel('time (fs)')
            
            plt.legend()
            if keep_eyield_vs_time == True:
                plt.savefig(f'IMAGE_yield_evolution_{variable_name}={var_i}.png')
            plt.close()

            # Anti - Normalization of the data
            if normalization == True:
                eyield = my_anti_renormalization(eyield, eyield_transform)
                t = my_anti_renormalization(t, t_transform)
                deyield_dt = my_anti_renormalization(deyield_dt, t_transform/eyield_transform)
                if t_eyield_end != None:
                    t_eyield_end = single_variable_renorm_to_antirenorm(t_eyield_end, t_transform)
        
            if t_eyield_end != None:
                data_eyield_end.append((var_i, t_eyield_end))
                #print(f'{prfx}For {variable_name} = {var_i}{variable_unit} : Cascade end is at {t_eyield_end} fs')
            else:
                print(f'{prfx}For {variable_name} = {var_i}{variable_unit} : Yield has not converged.')

            # Prepare for next step
            eyield = []
            t = []
            deyield_dt = []
    
    # Saving the data in a file
    data_eyield_end.sort(key=lambda x: x[0])
    with open(f'eyield_saturation_time_vs_{variable_name}.dat', 'w') as file:
        for i in range(len(data_eyield_end)):
            file.write(f'{data_eyield_end[i][0]}      {data_eyield_end[i][1]}\n')
    
    return data_eyield_end

# Math tools
def my_derivative(y=list,x=list):
    dy_dx = []
    Number_of_points = len(y)
        
    for i in range(Number_of_points):
    # general
        if i != 0 and i != (Number_of_points - 1): 
            dy_dx.append((y[i+1] - y[i-1])/(x[i+1]-x[i-1]))
    # first case
        if i == 0:
            dy_dx.append((y[i+1] - y[i])/(x[i+1]-x[i]))
        # last case
        if i == (Number_of_points-1): 
            dy_dx.append((y[i] - y[i-1])/(x[i]-x[i-1]))

    return dy_dx

def my_renormalization(x):
    length = len(x)
    x_renorm = x[-1] - x[0]
    for i in range(length):
        x[i] = x[i]/x_renorm
    return x, x_renorm

def my_anti_renormalization(x, x_renorm):
    length = len(x)
    for i in range(length):
        x[i] = x[i]*x_renorm
    return x

def single_variable_renorm_to_antirenorm(x_j, x_renorm):
    return x_j*x_renorm

# Settings printing:
def print_EYIELD_settings():
    print(f'{prfx}Analysis Settings')
    print(f'{prfx}Variable parameter: {variable_name}({variable_unit})')
    print(f'{prfx}Specific time analysis for: {time_to_read} fs')
    print(f'{prfx}Noise Filtering: {noise_filtering}')
    print(f'{prfx}Plot with Noise error bars (for yield): {plot_with_noise_error}')
    print(f'{prfx}Treshold percentage: {f*100}%')
    print(f'{prfx}Normalization: {normalization}')
    print(f'{prfx}Time reading error: {err_time} fs')
    print(f'{prfx}')
    print(f'{prfx}Production Settings')
    print(f'{prfx}Plot Electron Spectrum: {make_plot_e_spectrum}')
    print(f'{prfx}Plot Electron Yield: {make_plot_e_yield}')
    print(f'{prfx}Plot Electron Cascade Evolution End: {make_data_and_plot_end_of_cascades}')
    print(f'{prfx}Plot Electron Yield Evolution End: {make_data_and_plot_end_of_eyield}')
    print(f'{prfx}Make GIFs: {make_GIFs}')
    print(f'{prfx}')
    print(f'{prfx}Image Saving Settings')
    print(f'{prfx}Save Image >> {keep_e_spectrum_image} >> Electron Spectrum at {time_to_read} fs')
    print(f'{prfx}Save Image >> {keep_e_yield_image} >> Electron Yield at {time_to_read} fs')
    print(f'{prfx}Save Image >> {keep_cascade_time_vs_variable} >> Electron Cascade Evolution: Cascade Time VS {variable_name}')
    print(f'{prfx}Save Image >> {keep_electrons_vs_time} >> Electron Cascade Evolution: Electrons & Derivative VS Time')
    print(f'{prfx}Save Image >> {keep_eyield_saturation_vs_variable} >> Electron Yield Evolution: Saturation Time VS {variable_name}')
    print(f'{prfx}Save Image >> {keep_eyield_vs_time} >> Electron Yield Evolution: Yield & Derivative VS Time')
    print(f'{prfx}Save Image >> {keep_GIF_generated_images} >> Images to generate GIFs')

def print_subfolders_and_data(logical, directory_data_folders, number):
    print(f'{prfx}')
    if logical == True:
        if number > 1:
            print(f'{prfx}{directory_data_folders} contains {number} folders and they are:')
        else:
            print(f'{prfx}{directory_data_folders} contains {number} folder called:')
        counter = 1
        for subfolder in os.listdir(directory_data_folders):
            subfolder_path = os.path.join(directory_data_folders, subfolder)
            if os.path.isdir(subfolder_path):
                last_time = get_max_time(subfolder_path,file_e_spectrum)
                print(f'{prfx}{counter}) {subfolder}with d: {float(read_INPUT_DATA_material_z2(subfolder_path))}A, theta: {float(read_INPUT_DATA_theta(subfolder_path))}, barrier: {read_INPUT_DATA_barrier(subfolder_path)}, Time_sim: {last_time}fs')
                #print(f'{prfx}---> with d: {float(read_INPUT_DATA_material_z2(subfolder_path))}A, theta: {float(read_INPUT_DATA_theta(subfolder_path))}, barrier: {read_INPUT_DATA_barrier(subfolder_path)}, Time_sim: {last_time}fs')
            counter = counter + 1
    print(f'{prfx}')
    return round(last_time)       

# Program styling
def print_start_program(logical):
    if logical == True:
        print('************* E-YIELD: Data Analysis for TREKIS-4 ***********')
        print('')
        return 'EY >>> '
    else:
        return ''

def print_end_program(logical):
    if logical == True:
        print('')
        print('********************** End of E-YIELD ***********************')

# ************************************************************************
#                          START OF E-YIELD
# ************************************************************************

# Visual Output
prfx = print_start_program(Styling)

# Assigning settings for the variable to study
variable_name, variable_label, variable_unit = variable_settings(chosen_variable)

# Print the settings of EYIELD
print_EYIELD_settings()

# Checking if TREKIS outputs are present 
check_if_directory_exists(directory_data_folders)

# Counting the number of folders inside that directory
num_of_folders = count_the_subfolders(directory_data_folders)

# Get the time simulation & Printout the subfolders and check all their input data (significant variables)
time_simulation = print_subfolders_and_data(print_folders_and_significant_data, directory_data_folders, num_of_folders)


# ************************************************************************
#                          MAIN PROGRAM
# ************************************************************************

# Extract the The data and calculate the Electron Yield 
data = extract_data_from_spectrums() 

# Sort the variable in increasing order - necessary for the sorting algorithms
data.sort(key=lambda x: x[0]) 

# ************************************************************************
#                        PLOTTING AND GIF MAKING
# ************************************************************************

plot_electron_spectrum(data, time_to_read)

plot_electron_yield(data, time_to_read)

make_GIF_electron_spectrum_evolution(data)

make_GIF_electron_yield_evolution(data)

# ************************************************************************
#                          SATURATION ANALYSIS
# ************************************************************************

# Get the end time of cascade for each simulation
data_cascade_end = data_and_plot_end_of_cascades_vs_variable() # lists = ( (var1, time_cascade_end1), ...)
# UPDATE THIS: save, the data in a file, the images in a certain folder and add it as an option

# Sort the variable in increasing order - necessary for the sorting algorithms
data_cascade_end.sort(key=lambda x: x[0]) 

# Get the end time of electron yield saturation
data_e_yield_end = data_and_plot_end_of_eyield_vs_variable(data) # lists = ( (var1, time_eyield_end1), ...)
# UPDATE THIS: save, the data in a file, the images in a certain folder and add it as an option

# ************************************************************************
#                         PLAYGROUND E-YIELD
# ************************************************************************
x1, y1 = [], []
x2, y2 = [], []

material = 'Si'

for i in range(len(data_cascade_end)):
    x1.append(data_cascade_end[i][0])
    y1.append(data_cascade_end[i][1])

for i in range(len(data_e_yield_end)):
    x2.append(data_e_yield_end[i][0])
    y2.append(data_e_yield_end[i][1])

# Yield vs Electron Saturation 
plt.plot(x1, y1, marker='o', linestyle = '-', label = 'Saturation of $e^{-}$')
plt.plot(x2, y2, marker = 'o', linestyle = '-', label = 'Saturation of $e^{-}$ Yield')
plt.title(f'Cascade and Yield End Time VS {variable_name} in {material}')
plt.xlim(0)
plt.ylim(0)
plt.legend()
plt.xlabel(f'{variable_name} ({variable_unit})')
plt.ylabel(f'Cascade end time (fs)')
plt.savefig(f'IMAGE_cascade_end_VS_{variable_name}.png')
#plt.show()
plt.close()

# ************************************************************************
#                          Playground 2.0
# ************************************************************************
#                           The electron range
# ************************************************************************

data_e_density = extract_data_from_electron_density() # variable, time, coordinates, densities

# Finding the electron range

data_e_density.sort(key=lambda x: x[0]) #sort by  variable
data_e_density.sort(key=lambda x: x[1]) #sort by time

electron_energy = []
electron_range = []
electron_average_range = []

#print(data_e_density[-1])

#exit()

for i in range(len(data_e_density) - num_of_folders, len(data_e_density) ):

    #print(data_e_density[i])
    
    irradiation_energy = data_e_density[i][0]
    list_coordinates = data_e_density[i][2]
    list_densities = data_e_density[i][3]

    max_electron_coordinate = 0

    for j in range(len(list_coordinates)):
        if list_densities[j] !=0 and list_coordinates[j] > max_electron_coordinate:
            max_electron_coordinate = list_coordinates[j]
    
    electron_energy.append(irradiation_energy) # the var1 - which in this case would be the irradiation energy always
    electron_range.append(max_electron_coordinate/10) # division by 10 to pass to nanometers

    # Calculate the weighted average range
    list_weighted_z_coordinate = []
    for j in range(len(list_coordinates)):
        list_weighted_z_coordinate.append(list_coordinates[j]*list_densities[j])

    integral_1 = scipy.integrate.simpson(list_coordinates, list_weighted_z_coordinate)
    integral_2 = scipy.integrate.simpson(list_coordinates, list_densities)

    averaged_coordinate = integral_1/integral_2

    if integral_2 != 0:
        electron_average_range.append(averaged_coordinate/10) # the 10 is here to switch from A to nm
    else:
        electron_average_range.append(0)


    print(f'{irradiation_energy} {max_electron_coordinate}')
    print(f'averaged result is {averaged_coordinate}')

plt.plot(electron_energy, electron_range, marker = 'o', linestyle = '-', color = 'green', label = 'max')
plt.plot(electron_energy, electron_average_range, marker = 'o', linestyle = '-', color = 'purple', label = 'average')
plt.xlabel('Electron energy(eV)')
plt.ylabel('Electron range(nm)')
plt.xlim(1e1,1e5)
plt.ylim(1e0, 1e4)
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
plt.close()

# ************************************************************************
#                          END OF E-YIELD
# ************************************************************************

print_end_program(Styling)
