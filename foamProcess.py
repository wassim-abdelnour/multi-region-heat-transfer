import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import figure

def print_Process_log(Process_log, sim_log):
    if(Process_log == '\n'):
        print()
        sim_log.write(Process_log)
    else:
        print(Process_log)
        sim_log.write(Process_log + '\n')

def compute_head_losses(simulation_name, number_of_iterations, Px, Py, Pz, sim_log):
    os.system('mpirun -np ' + str(Px*Py*Pz) + ' postProcess -case ' + simulation_name + ' -latestTime -parallel -region fluid -func \"patchAverage(fields=(p), patch=Inlet)\" > ' + simulation_name + '/Head_Losses.log')
    x = np.genfromtxt(simulation_name + '/postProcessing/fluid/patchAverage(fields=(p),patch=Inlet)/' + str(number_of_iterations) + '/surfaceFieldValue.dat', delimiter='\t')
    head_losses = x[1]
    print_Process_log('head losses between Inlet and Outlet are : ' + "{:e}".format(head_losses) + ' Pa', sim_log)
    print_Process_log('\n', sim_log)
    return head_losses

def compute_Yplus(simulation_name, number_of_iterations, Px, Py, Pz, sim_log):
    os.system('mpirun -np ' + str(Px*Py*Pz) + ' chtMultiRegionFoam -case ' + simulation_name + ' -postProcess -parallel -region fluid -func yPlus > ' + simulation_name + '/Yplus.log')

    x = np.genfromtxt(simulation_name + '/postProcessing/fluid/yPlus/' + str(number_of_iterations) + '/yPlus.dat', delimiter='\t')
    yPlus_avg = x[4] #ligne #4colonne
    yPlus_max = x[3]

    print_Process_log('average Y+ is : ' + "{:e}".format(yPlus_avg), sim_log)
    print_Process_log('maximum Y+ is : ' + "{:e}".format(yPlus_max), sim_log)
    print_Process_log('\n', sim_log)

    return yPlus_avg, yPlus_max

def compute_thermal_resistance(simulation_name, number_of_iterations, Q, Px, Py, Pz, sim_log):
    os.system('mpirun -np ' + str(Px*Py*Pz) + ' postProcess -case ' + simulation_name + ' -latestTime -region fluid -parallel -func \"patchAverage(fields=(T), patch=Inlet)\" > ' + simulation_name + '/Tin_avg.log')
    x = np.genfromtxt(simulation_name + '/postProcessing/fluid/patchAverage(fields=(T),patch=Inlet)/' + str(number_of_iterations) + '/surfaceFieldValue.dat', delimiter='\t')
    Tin_avg = x[1]

#mpirun -np 9 postProcess -case 01_debit_0.008700 -latestTime -region fluid -parallel -func "patchAverage(fields=(T), patch=Outlet)"
    os.system('mpirun -np ' + str(Px*Py*Pz) + ' postProcess -case ' + simulation_name + ' -latestTime -region fluid -parallel -func \"patchAverage(fields=(T), patch=Outlet)\" > ' + simulation_name + '/Tout_avg.log')
    x = np.genfromtxt(simulation_name + '/postProcessing/fluid/patchAverage(fields=(T),patch=Outlet)/' + str(number_of_iterations) + '/surfaceFieldValue.dat', delimiter='\t')
    Tout_avg = x[1]

    os.system('mpirun -np ' + str(Px*Py*Pz) + ' postProcess -case ' + simulation_name + ' -latestTime -region solid -parallel -func \"patchAverage(fields=(T), patch=Heated_Face)\" > ' + simulation_name + '/Th_avg.log')
    x = np.genfromtxt(simulation_name + '/postProcessing/solid/patchAverage(fields=(T),patch=Heated_Face)/' + str(number_of_iterations) + '/surfaceFieldValue.dat', delimiter='\t')
    Th_avg = x[1]

    #thermal_resistance = [ Theat - (Tin+Tout)/2 ] / Q
    thermal_resistance = ( Th_avg - Tin_avg ) / Q


    print_Process_log('thermal resistance is : ' + "{:e}".format(thermal_resistance) + ' K/W', sim_log)
    print_Process_log('\n', sim_log)

    return thermal_resistance

def plot_convergence_curve(simulation_name, number_of_iterations, pack_of_iterations, Q, rho, Cp_f, debit_massique, Tin, sim_log):
    plot_interface_heat_flux_convergence(simulation_name, number_of_iterations, pack_of_iterations, Q, sim_log)
    plot_outlet_heat_flux_convergence(simulation_name, number_of_iterations, pack_of_iterations, Q, rho, Cp_f, debit_massique, Tin, sim_log)
    plot_outlet_mass_flow_rate_convergence(simulation_name, number_of_iterations, pack_of_iterations, debit_massique, sim_log)

def plot_interface_heat_flux_convergence(simulation_name, number_of_iterations, pack_of_iterations, Q, sim_log):
    print_Process_log('plotting convergence curve 1 ...', sim_log)
    main_list = np.genfromtxt(simulation_name + '/postProcessing/fluid/interfaceHeatFlux/' + '0' + '/wallHeatFlux.dat', delimiter='\t')
    for iter in range(pack_of_iterations, number_of_iterations, pack_of_iterations):
        auxiliary_list = np.genfromtxt(simulation_name + '/postProcessing/fluid/interfaceHeatFlux/' + str(iter) + '/wallHeatFlux.dat', delimiter='\t')
        main_list = np.concatenate((main_list, auxiliary_list))
    main_list = np.delete(main_list, 1, 1)
    plt.figure(figsize =(6, 4), dpi = 1500)
    plt.title("Convergence Plot 1")
    plt.xlabel("Number of Iterations")
    plt.ylabel("Heat Flux Through Interface (W)")
    plt.grid(linestyle = '--', linewidth = 0.1)
    plt.plot(main_list[:, 0], main_list[:,3], 'b-', label = 'At Interface')
    plt.plot(main_list[:, 0], [Q for x in range(len(main_list[:, 0]))], 'r-', label = 'At Source')
    plt.legend(loc=4)
    plt.savefig(simulation_name + '/interface_heat_flux_convergence_plot.png')
    print_Process_log('convergence plot saved as png', sim_log)
    print_Process_log('\n', sim_log)

def plot_outlet_heat_flux_convergence(simulation_name, number_of_iterations, pack_of_iterations, Q, rho, Cp_f, debit_massique, Tin, sim_log):
    print_Process_log('plotting convergence curve 2 ...', sim_log)
    main_list = np.genfromtxt(simulation_name + '/postProcessing/fluid/heatfluxoutlet/' + '0' + '/fieldValueDelta.dat', delimiter='\t')
    for iter in range(pack_of_iterations, number_of_iterations, pack_of_iterations):
        auxiliary_list = np.genfromtxt(simulation_name + '/postProcessing/fluid/heatfluxoutlet/' + str(iter) + '/fieldValueDelta.dat', delimiter='\t')
        main_list = np.concatenate((main_list, auxiliary_list))
    main_list[:,1][main_list[:,1]<0] = 0
    plt.figure(figsize =(6, 4), dpi = 1500)
    plt.title("Convergence Plot 2")
    plt.xlabel("Number of Iterations")
    plt.ylabel("Heat Flux Through Outlet (W)")
    plt.grid(linestyle = '--', linewidth = 0.1)
    plt.plot(main_list[:, 0], main_list[:,1], 'b-', label = 'At Outlet')
    plt.plot(main_list[:, 0], [Q for x in range(len(main_list[:, 0]))], 'r-', label = 'At Source')
    plt.legend(loc=4)
    plt.savefig(simulation_name + '/outlet_heat_flux_convergence_plot.png')
    print_Process_log('convergence plot saved as png', sim_log)
    print_Process_log('\n', sim_log)

def plot_outlet_mass_flow_rate_convergence(simulation_name, number_of_iterations, pack_of_iterations, debit_massique, sim_log):
    print_Process_log('plotting convergence curve 3 ...', sim_log)
    main_list = np.genfromtxt(simulation_name + '/postProcessing/fluid/outletflowrate/' + '0' + '/surfaceFieldValue.dat', delimiter='\t')
    for iter in range(pack_of_iterations, number_of_iterations, pack_of_iterations):
        auxiliary_list = np.genfromtxt(simulation_name + '/postProcessing/fluid/outletflowrate/' + str(iter) + '/surfaceFieldValue.dat', delimiter='\t')
        main_list = np.concatenate((main_list, auxiliary_list))
    plt.figure(figsize =(6, 4), dpi = 1500)
    plt.title("Convergence Plot 3")
    plt.xlabel("Number of Iterations")
    plt.ylabel("Mass Flow Rate Through Outlet (kg/s)")
    plt.grid(linestyle = '--', linewidth = 0.1)
    plt.plot(main_list[:, 0], main_list[:,1], 'b-', label = 'At Outlet')
    plt.plot(main_list[:, 0], [debit_massique for x in range(len(main_list[:, 0]))], 'r-', label = 'At Inlet')
    plt.legend(loc=4)
    plt.savefig(simulation_name + '/outlet_mass_flow_rate_convergence_plot.png')
    print_Process_log('convergence plot saved as png', sim_log)
    print_Process_log('\n', sim_log)

def plot_residuals(simulation_name, number_of_iterations, pack_of_iterations, sim_log):
    iterations_list = [x for x in range(1, number_of_iterations+1, 1)]
    print_Process_log('plotting residuals ...', sim_log)
    residuals_Ux = np.genfromtxt(simulation_name + '/logs/UxFinalRes_0', delimiter='\t')
    print_Process_log('final residual of Ux plotted', sim_log)
    residuals_Uy = np.genfromtxt(simulation_name + '/logs/UyFinalRes_0', delimiter='\t')
    print_Process_log('final residual of Uy plotted', sim_log)
    residuals_Uz = np.genfromtxt(simulation_name + '/logs/UzFinalRes_0', delimiter='\t')
    print_Process_log('final residual of Uz plotted', sim_log)
    residuals_hf = np.genfromtxt(simulation_name + '/logs/hFinalRes_0', delimiter='\t')
    print_Process_log('final residual of hf plotted', sim_log)
    residuals_hs = np.genfromtxt(simulation_name + '/logs/eFinalRes_0', delimiter='\t')
    print_Process_log('final residual of es plotted', sim_log)
    plt.clf()
    plt.figure(figsize =(6, 4), dpi = 1500)
    plt.title("Final Resiudals")
    plt.xlabel("Number of Iterations")
    plt.ylabel("Residual")
    plt.grid(linestyle = '--', linewidth = 0.1)
    plt.plot(iterations_list, residuals_Ux[:,1], 'b-', label = 'Ux')
    plt.plot(iterations_list, residuals_Uy[:,1], 'g-', label = 'Uy')
    plt.plot(iterations_list, residuals_Uz[:,1], 'y-', label = 'Uz')
    plt.plot(iterations_list, residuals_hf[:,1], 'k-', label = 'hf')
    plt.plot(iterations_list, residuals_hs[:,1], 'm-', label = 'hs')
    plt.legend(loc=1)
    plt.savefig(simulation_name + '/residuals_plot.png')
    for iter in range(0, number_of_iterations, pack_of_iterations):
        plt.clf()
    	plt.figure(figsize =(6, 4), dpi = 1500)
    	plt.title("Final Resiudals")
    	plt.xlabel("Number of Iterations")
    	plt.ylabel("Residual")
    	plt.grid(linestyle = '--', linewidth = 0.1)
        plt.plot(iterations_list[iter:iter+pack_of_iterations], residuals_Ux[iter:iter+pack_of_iterations,1], 'b-', label = 'Ux')
        plt.plot(iterations_list[iter:iter+pack_of_iterations], residuals_Uy[iter:iter+pack_of_iterations,1], 'g-', label = 'Uy')
        plt.plot(iterations_list[iter:iter+pack_of_iterations], residuals_Uz[iter:iter+pack_of_iterations,1], 'y-', label = 'Uz')
        plt.plot(iterations_list[iter:iter+pack_of_iterations], residuals_hf[iter:iter+pack_of_iterations,1], 'k-', label = 'hf')
        plt.plot(iterations_list[iter:iter+pack_of_iterations], residuals_hs[iter:iter+pack_of_iterations,1], 'm-', label = 'hs')
        plt.legend(loc=1)
        plt.savefig(simulation_name + '/residuals_plot_' + str(iter) + '-' + str(iter+pack_of_iterations) + '.png')              
    print_Process_log('residuals plot saved as png', sim_log)
    print_Process_log('\n', sim_log)

def compute_simulation_time(simulation_name, number_of_iterations, pack_of_iterations, sim_log):
    time_list=[]
    x = np.genfromtxt(simulation_name + '/logs/executionTime_0', delimiter='\t')
    for i in range(1, number_of_iterations+1, 1):
        if(i%pack_of_iterations == 0):
            time_list.append(x[i-1][1])
    simulation_time = sum(time_list)
    print_Process_log('total simulation time is : ' + "{0:.2f}".format(simulation_time) + ' s', sim_log)
    print_Process_log('\n', sim_log)
