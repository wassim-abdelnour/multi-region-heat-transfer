import math
import os
import foamSc
import numpy as np

def print_Ctrl_log(Ctrl_log, sim_log):
    if(Ctrl_log == '\n'):
        print()
        sim_log.write(Ctrl_log)
    else:
        print(Ctrl_log)
        sim_log.write(Ctrl_log + '\n')

def get_simulation_name(i, debit_massique):
    if i < 10 :
        prefix = '0'
    else:
        prefix = ''
    simulation_name = prefix + str(i) + '_debit_' + "{0:.6f}".format(debit_massique)
    os.system('mkdir ' + simulation_name)

    return simulation_name

def print_simulation_recap(plate_width, plate_height, plate_thickness, longueur, largeur, spacing, perc_of_r, nombre_de_passes, debit_massique, rho, mu, Tin, Cp_f, kf, Q, rho_s, Cp_s, ks, Px, Py, Pz, sim_log):

    print_Ctrl_log('''
-------------------------- Temisth --------------------------
***Technology Leader for High Performances Thermal systems***
-------------------------------------------------------------

 Cold Plates Thermo-Physical Caracterization using openFoam 

-------------------------------------------------------------

 version 3.0 - 7 Fev 2023

-------------------------------------------------------------

''', sim_log)

    ####################Print Simulation Recap####################
    print_Ctrl_log('plate dimensions : ' + str(plate_width) + ' x ' + str(plate_height) + ' x ' + str(plate_thickness) + ' m', sim_log)
    print_Ctrl_log('duct dimensions : ' + str(longueur) + ' x ' + str(largeur) + ' m', sim_log)
    print_Ctrl_log('space between duct and plate walls : ' + str(spacing) + ' m', sim_log)
    print_Ctrl_log('number of passes : ' + str(nombre_de_passes) + ', with corner radius = duct_width*' + str(perc_of_r) + '%', sim_log)
    print_Ctrl_log('\n', sim_log)
    print_Ctrl_log('mass flow rate : ' + "{0:.6f}".format(debit_massique) + ' kg/s', sim_log)
    print_Ctrl_log('fluid density : ' + str(rho) + ' kg/m3', sim_log)
    print_Ctrl_log('fluid dynamic viscosity : ' + str(mu) + ' Pa.s', sim_log)
    print_Ctrl_log('fluid temperature at inlet : ' + str(Tin) + ' K', sim_log)
    print_Ctrl_log('fluid specific heat : ' + str(Cp_f) + ' J/kg.K', sim_log)
    print_Ctrl_log('fluid conductivity : ' + str(kf) + ' W/m.K', sim_log)
    print_Ctrl_log('\n', sim_log)
    print_Ctrl_log('heat flux on heated wall : ' + str(Q) + ' W', sim_log)
    print_Ctrl_log('solid density : ' + str(rho_s) + ' kg/m3', sim_log)
    print_Ctrl_log('solid specific heat : ' + str(Cp_s) + ' J/kg.K', sim_log)
    print_Ctrl_log('solid conductivity : ' + str(ks) + ' W/m.K', sim_log)
    print_Ctrl_log('\n', sim_log)
    print_Ctrl_log('Decomposition for Parallelization (width x height x thickness) : ' + str(Px) + ' x ' + str(Py) + ' x ' + str(Pz), sim_log)
    print_Ctrl_log('\n', sim_log)
    ##############################################################

def compute_turbulence_parameters(longueur, largeur, debit_massique, rho, mu, sim_log):
    #inlet area = Pi*r2
    S = longueur * largeur
    P = 2 * (longueur + largeur)
    Hd = 4*S / P

    #debit_massique = rho * V * S
    Vmoyenne = debit_massique / (rho*S)

    #free-stream velocity = vitesse maximale d'un profile paraboloique au meme débit
    Vinf = 2*Vmoyenne #à revoir, puisque le tube n'est plus symétrique. Mais bonne approximation pour l'instent

    ################Compute Turbulence Parameters#################
    #The caracteristic length for Re calculation is the hydraulic diameter
    #Hydraulic Diameter = 4*Surface_Area/Perimeter
    #Hydraulic Number of a Rectangular cross section is 4*(L*l) / 2*(L+l)
    Re = rho * Vmoyenne * Hd / mu
    print_Ctrl_log('Re = ' + str(int(Re)), sim_log)
    print_Ctrl_log('\n', sim_log)

    #Calculation of k (kinetic Energy)
    I = 0.16*pow(Re, -1/8) #Turbulence Intensity
    k = 1.5 * Vinf*Vinf * I*I #kappa

    #Calculaction of epsilon (turbulent dissipation rate)
    l = 0.07 * Hd #mixing length or turbulent length scale
    eps = (math.pow(0.09,0.75)*math.pow(k,1.5)) / (l)

    #Calculaction of nut (turbulent viscosity)
    nut = 0.09 * k*k/eps

    #Calculaction of omega (specific turbulent dissipation rate)
    omega = math.pow(0.09, -0.25) * math.sqrt(k)/l

    return k, eps, nut, omega, l

def run_simulation(simulation_name, pack_of_iterations, Q, debit_massique, Cp_f, convergence_criterion, Px, Py, Pz, sim_log):
    number_of_iterations = pack_of_iterations
    print_Ctrl_log('Decomposing Mesh into Regions for Parallellization ...', sim_log)
    os.system('decomposePar -case ' + simulation_name + ' -allRegions > ' + simulation_name + '/decomposePar.log')
    print_Ctrl_log('Decomposing Completed', sim_log)
    print_Ctrl_log('\n', sim_log)
    print_Ctrl_log('Simulation is running ...', sim_log)

    converged = False
    while(converged == False):
        os.system('mpirun -np ' + str(Px*Py*Pz) + ' chtMultiRegionFoam -case ' + simulation_name + ' -parallel > ' + simulation_name + '/solverLogs/solver_' + str(number_of_iterations-pack_of_iterations) + '.log')

        heat_flux_interface_array = np.genfromtxt(simulation_name + '/postProcessing/fluid/interfaceHeatFlux/' + str(number_of_iterations-pack_of_iterations) + '/wallHeatFlux.dat', delimiter='\t')
        heat_flux_interface_array = np.delete(heat_flux_interface_array, 1, 1) #to clean nan column
        mean_heat_flux_interface = np.mean(heat_flux_interface_array, axis=0)[3]
        error_heat_flux_interface = abs(mean_heat_flux_interface - Q )*100 / Q
        print_Ctrl_log('error on thermal transfer through interface for iterations ' + str(number_of_iterations-pack_of_iterations) + '-' + str(number_of_iterations) + ' is : ' + "{0:.2f}".format(error_heat_flux_interface) +' %', sim_log)

        heat_flux_outlet_array = np.genfromtxt(simulation_name + '/postProcessing/fluid/heatfluxoutlet/' + str(number_of_iterations-pack_of_iterations) + '/fieldValueDelta.dat', delimiter='\t')
        heat_flux_outlet_array = np.delete(heat_flux_outlet_array, 0, 0) #to not compute zro values in the error
        heat_flux_outlet = np.mean(heat_flux_outlet_array, axis=0)[1]
        error_heat_flux_outlet = abs(heat_flux_outlet - Q )*100 / Q #il faut le multiplier par RhoCp et le soustraire par MpCpT entree puis le comparer avec Q
        print_Ctrl_log('error on heat flux through outlet for iterations ' + str(number_of_iterations-pack_of_iterations) + '-' + str(number_of_iterations) + ' is : ' + "{0:.2f}".format(error_heat_flux_outlet) +' %', sim_log)

        mass_flow_outlet_array = np.genfromtxt(simulation_name + '/postProcessing/fluid/outletflowrate/' + str(number_of_iterations-pack_of_iterations) + '/surfaceFieldValue.dat', delimiter='\t')
        mass_flow_outlet = np.mean(mass_flow_outlet_array, axis=0)[1]
        error_mass_flow_outlet = abs(mass_flow_outlet - debit_massique )*100 / debit_massique
        print_Ctrl_log('error on mass flow rate through outlet for iterations ' + str(number_of_iterations-pack_of_iterations) + '-' + str(number_of_iterations) + ' is : ' + "{0:.2f}".format(error_mass_flow_outlet) +' %', sim_log)

        print_Ctrl_log('\n', sim_log)
        if(error_heat_flux_interface < convergence_criterion and error_heat_flux_outlet < convergence_criterion and error_mass_flow_outlet < convergence_criterion):
            converged = True
        else:
            for procNum in range(0, Px*Py*Pz, 1):
                os.system('rm -r ' + simulation_name + '/processor' + str(procNum) + '/' + str(number_of_iterations-pack_of_iterations))
            number_of_iterations += pack_of_iterations
            foamSc.write_Control_Dict(simulation_name, number_of_iterations, pack_of_iterations, Cp_f)
    print_Ctrl_log('Simulation Ended', sim_log)
    print_Ctrl_log('\n', sim_log)

    return number_of_iterations

def reconstruct_simulation(simulation_name, sim_log):
    print_Ctrl_log('recomposing mesh partitions ...', sim_log)
    os.system('reconstructPar -case ' + simulation_name + ' -allRegions > ' + simulation_name + '/reconstruction.log')
    print_Ctrl_log('recomposition completed', sim_log)
    os.system('rm -r ' + simulation_name +'/processor*')
    print_Ctrl_log('partions cleared', sim_log)
    print_Ctrl_log('\n', sim_log)

def Extract_Solver_Details(simulation_name, number_of_iterations, pack_of_iterations, sim_log):
    solver_log = open('solver.log', "w")
    print_Ctrl_log('merging solvers log files ...', sim_log)
    for x in range(0, number_of_iterations, pack_of_iterations):
        sub_solver_log = open(simulation_name + '/solverLogs/solver_' + str(x) +'.log', "r")
        sub_solver_data = sub_solver_log.read()
        solver_log.write(sub_solver_data)
        sub_solver_log.close()
    solver_log.close()
    print_Ctrl_log('solvers log files merged', sim_log)
    print_Ctrl_log('\n', sim_log)
    os.system('rm -r ' + simulation_name + '/solverLogs')
    print_Ctrl_log('extracting solver details ...', sim_log)
    os.system('foamLog solver.log > ' + simulation_name + '/solverLogs.log')
    os.system('rm logs/clockTime_0')
    os.system('rm logs/contCumulative_0')
    os.system('rm logs/contGlobal_0')
    os.system('rm logs/contLocal_0')
    os.system('rm logs/fluid_0')
    os.system('rm logs/fluidFinalRes_0')
    os.system('rm logs/fluidIters_0')
    os.system('rm logs/h_0')
    os.system('rm logs/e_0')
    os.system('rm logs/hIters_0')
    os.system('rm logs/eIters_0')
    os.system('rm logs/k_0')
    os.system('rm logs/kIters_0')
    os.system('rm logs/omega_0')
    os.system('rm logs/p_rgh_0')
    os.system('rm logs/p_rghIters_0')
    os.system('rm logs/Separator_0')
    os.system('rm logs/solid_0')
    os.system('rm logs/solidFinalRes_0')
    os.system('rm logs/solidIters_0')
    os.system('rm logs/Time_0')
    os.system('rm logs/Ux_0')
    os.system('rm logs/UxIters_0')
    os.system('rm logs/Uy_0')
    os.system('rm logs/UyIters_0')
    os.system('rm logs/Uz_0')
    os.system('rm logs/UzIters_0')
    os.system('rm logs/CourantMax_0')
    os.system('rm logs/CourantMax_1')
    os.system('rm logs/CourantMean_0')
    os.system('rm logs/CourantMean_1')
    os.system('mv logs ' + simulation_name)
    os.system('mv solver.log ' + simulation_name)
    print_Ctrl_log('solver details extracted', sim_log)
    print_Ctrl_log('\n', sim_log)
