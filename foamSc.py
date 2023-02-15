import os

def create_openFoam_simulation_setup(simulation_name, debit_massique, rho, mu, Cp_f, kf, Q, Tin, rho_s, Cp_s, ks, Px, Py, Pz, k, eps, nut, omega, l, number_of_iterations, pack_of_iterations):
#   os.system('mkdir ' + simulation_name)
    os.system('mkdir ' + simulation_name + '/0')
    os.system('mkdir ' + simulation_name + '/0/fluid')
    os.system('mkdir ' + simulation_name + '/0/solid')
    os.system('mkdir ' + simulation_name + '/constant')
    os.system('mkdir ' + simulation_name + '/constant/fluid')
    os.system('mkdir ' + simulation_name + '/constant/solid')
    os.system('mkdir ' + simulation_name + '/system')
    os.system('mkdir ' + simulation_name + '/system/fluid')
    os.system('mkdir ' + simulation_name + '/system/solid')
    os.system('mkdir ' + simulation_name + '/solverLogs')

    write_gravity_file(simulation_name)
    write_region_properties(simulation_name)
    write_fluid_radiation_properties(simulation_name)
    write_turbulence_properties(simulation_name)
    write_fluid_thermophysical_properties(simulation_name, rho, Cp_f, kf, mu)
    write_solid_radiation_properties(simulation_name)
    write_solid_thermophysical_properties(simulation_name, rho_s, Cp_s, ks)
    write_fvSolution(simulation_name)
    write_Parallel_Decomposition(simulation_name, Px, Py, Pz)
    write_Control_Dict(simulation_name, number_of_iterations, pack_of_iterations, Cp_f)
    write_fluid_fvSchemes(simulation_name)
    write_fluid_fvSolution(simulation_name)
    write_solid_fvSchemes(simulation_name)
    write_solid_fvOptions(simulation_name)
    write_solid_fvSolution(simulation_name)
    write_solid_temperature(simulation_name, Tin, Q)
    write_fluid_velocity(simulation_name, debit_massique, rho)
    write_fluid_pressure(simulation_name)
    write_fluid_pressure_rgh(simulation_name)
    write_fluid_temperature(simulation_name, Tin)
    write_alphat(simulation_name)
    write_turbulent_intensity(simulation_name, k)
    write_turbulent_viscosity(simulation_name, nut)
    write_turbulent_dissipation_rate(simulation_name, omega, l)
    write_density(simulation_name, rho)
    
    os.system('touch ' + simulation_name + '/foam.foam')

header = '''
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                  |                   Temisth                      |
|    ||      T  ransferts    |    technology leader for high performances     |
|    ||      E  nergétiques  |                thermal systems                 |
|    ||      M  atériaux     |                                                |
| ___        I  nnovants     |    Web:      http://temisth.com/               |
|  |   |_    S  ystèmes      |                                                |
|  |   | |   Th ermiques     | This file was auto generated, do not modify it |
\*---------------------------------------------------------------------------*/'''

def write_gravity_file(simulation_name):
    f = open(simulation_name + "/constant/fluid/g", "w")
    gravity_file = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       uniformDimensionedVectorField;
    object      g;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -2 0 0 0 0];
value           (0 0 0);

// ************************************************************************* //
'''
    f.write(gravity_file)
    f.close()

def write_region_properties(simulation_name):
    f = open(simulation_name + "/constant/regionProperties", "w")
    region_properties =  header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      regionProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

regions
(
    fluid       (fluid)
    solid       (solid)
);

// ************************************************************************* //
'''
    f.write(region_properties)
    f.close()

def write_fluid_radiation_properties(simulation_name):
    f = open(simulation_name + "/constant/fluid/radiationProperties", "w")
    fluid_radiation_properties =  header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      radiationProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

radiation       off;

radiationModel  none;


// ************************************************************************* //
'''
    f.write(fluid_radiation_properties)
    f.close()

def write_turbulence_properties(simulation_name):
    f = open(simulation_name + "/constant/fluid/turbulenceProperties", "w")
    turbulence_properties =  header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType  RAS;

RAS
{
    RASModel        kOmegaSST;

    turbulence      on;

    printCoeffs     on;
}

// ************************************************************************* //
'''
    f.write(turbulence_properties)
    f.close()

def write_fluid_thermophysical_properties(simulation_name, rho, Cp_f, kf, mu):
    Pr = Cp_f * mu / kf
    f = open(simulation_name + "/constant/fluid/thermophysicalProperties", "w")
    fluid_thermophysical_properties = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}

mixture
{
    specie
    {
        molWeight       18;
    }
    equationOfState
    {
        rho    ''' + str(rho) + ''';
    }
    thermodynamics
    {
        Hf		0;
        Sf		0;
        Cp		''' + str(Cp_f) + ''';
    }
    transport
    {
        mu		''' + str(mu) + ''';
	Pr		''' + str(Pr) + ''';
    }
}

// ************************************************************************* //
'''
    f.write(fluid_thermophysical_properties)
    f.close()

def write_solid_radiation_properties(simulation_name):
    f = open(simulation_name + "/constant/solid/radiationProperties", "w")
    solid_radiation_properties = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      radiationProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

radiation       off;

radiationModel  none;


// ************************************************************************* //
'''
    f.write(solid_radiation_properties)
    f.close()

def write_solid_thermophysical_properties(simulation_name, rho_s, Cp_s, ks):
    f = open(simulation_name + "/constant/solid/thermophysicalProperties", "w")
    solid_thermophysical_properties = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

thermoType
{
    type            heSolidThermo;
    mixture         pureMixture;
    transport       constIsoSolid;
    thermo          eConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleInternalEnergy;
}

mixture
{
    specie
    {
        molWeight   26.98;
    }

    transport
    {
        kappa   '''+ str(ks) +''';
    }

    thermodynamics
    {
        Hf      0;
        Cv      '''+ str(Cp_s) +''';
    }

    equationOfState
    {
        rho     '''+ str(rho_s) +''';
    }
}

// ************************************************************************* //
'''
    f.write(solid_thermophysical_properties)
    f.close()

def write_fvSolution(simulation_name):
    f = open(simulation_name + "/system/fvSolution", "w")
    fvSolution = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

PIMPLE
{
    nOuterCorrectors 1;
}

// ************************************************************************* //

'''
    f.write(fvSolution)
    f.close()

def write_Parallel_Decomposition(simulation_name, Px, Py, Pz):
    nbre_de_coeurs = Px*Py*Pz #idéalement Pz <= 2
    f = open(simulation_name + "/system/decomposeParDict", "w")
    Parallel_Decomposition = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains  ''' + str(nbre_de_coeurs) + ''';

method          scotch;

coeffs
{
    n           (''' + str(Px) + ' ' + str(Py) + ' ' + str(Pz) + ''');
}

scotchCoeffs
{
}

// ************************************************************************* //
'''
    f.write(Parallel_Decomposition)
    f.close()

def write_Control_Dict(simulation_name, number_of_iterations, pack_of_iterations, Cp_f):
    writeInterval = pack_of_iterations/10 #100
    f = open(simulation_name + "/system/controlDict", "w")
    Control_Dict = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     chtMultiRegionFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         ''' + str(number_of_iterations) + ''';

deltaT          1;

writeControl    timeStep;

writeInterval   ''' + str(int(writeInterval)) +'''; //50

purgeWrite      1; //2

writeFormat     ascii;

writePrecision  7;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable yes;

maxCo           0.6;

// Maximum diffusion number
maxDi           10.0;

adjustTimeStep  no;

functions
{
    interfaceHeatFlux
    {
	// Mandatory entries (unmodifiable)
	type            wallHeatFlux;
	libs            ("libfieldFunctionObjects.so");

	// Optional entries (runtime modifiable)
	patches     (fluid_to_solid); // (wall1 "(wall2|wall3)");
	qr          qr;

	// Optional (inherited) entries
	writePrecision  8;
	writeToFile     true;
	useUserTime     true;
	region          fluid;
	enabled         true;
	log             true;
	timeStart       0;
//	timeEnd         5000;
	executeControl  timeStep;
	executeInterval ''' + str(int(writeInterval)) +''';
	writeControl    timeStep;
	writeInterval   ''' + str(int(writeInterval)) +''';
    }

    Co
    {
        // Mandatory entries (unmodifiable)
        type            CourantNo;
        libs            ("libfieldFunctionObjects.so");

        // Optional entries (runtime modifiable)
        rho             rho;

        // Optional (inherited) entries
        field           phi;
        result          Co;
        region          fluid;
        enabled         true;
        log             true;
        timeStart       0;
//      timeEnd         1000;
        executeControl  timeStep;
        executeInterval ''' + str(int(writeInterval)) +''';
        writeControl    timeStep;
        writeInterval   ''' + str(int(writeInterval)) +''';
    }

    PecletNo
    {
	// Mandatory entries (unmodifiable)
        type            PecletNo;
        libs            ("libfieldFunctionObjects.so");

	// Optional entries (runtime modifiable)
        rho             rho;

	// Optional (inherited) entries
        field           phi;
        result          Pe;
        region          fluid;
        enabled         true;
        log             true;
        timeStart       0;
//      timeEnd         1000;
        executeControl  timeStep;
        executeInterval ''' + str(int(writeInterval)) +''';
        writeControl    timeStep;
        writeInterval   ''' + str(int(writeInterval)) +''';
    }

    NormalVelocity
    {
        type            components;
        libs            ("libfieldFunctionObjects.so");

        field           U;

        region          fluid;
        enabled         true;
        log             true;
        timeStart       0;
        executeControl  timeStep;
        executeInterval ''' + str(int(writeInterval)) +''';
        writeControl    timeStep;
        writeInterval   ''' + str(int(writeInterval)) +''';
    }

   Cp
    {
        type            uniform;
        libs            ("libfieldFunctionObjects.so");

        fieldType       volScalarField;
        name            Cp_f;
        dimensions      [0 2 -2 -1 0 0 0];
        value           ''' + str(-Cp_f) +''';
        region		fluid;

        enabled         true;
        log             true;
        timeStart       0;
        executeControl  timeStep;
        executeInterval ''' + str(int(writeInterval)) +''';
        writeControl    timeStep;
        writeInterval   ''' + str(int(writeInterval)) +''';
    }

    heatfluxfields
    {
        type            multiply;
        libs            ("libfieldFunctionObjects.so");

        fields          (Uy T Cp_f rho);

        result          heatfluxfields;
        region          fluid;
        enabled         true;
        log             true;
        timeStart       0;
        executeControl  timeStep;
        executeInterval ''' + str(int(writeInterval)) +''';
        writeControl    timeStep;
        writeInterval   ''' + str(int(writeInterval)) +''';
    }

    heatfluxoutlet
    {
        type            fieldValueDelta;
        libs            ("libfieldFunctionObjects.so");

        writeControl    timeStep;
        writeInterval   ''' + str(int(writeInterval)) +''';

        operation       add;
        region fluid    fluid;
        writeFields     true;
        log             true;

        region1
        {
            type            surfaceFieldValue;
            libs            ("libfieldFunctionObjects.so");

            writeControl    timeStep;
            writeInterval   ''' + str(int(writeInterval)) +''';

            writeFields     false;
            operation       areaIntegrate;
            fields          (heatfluxfields);
            regionType      patch;
            name            Inlet;
            region          fluid;
        }

        region2
        {
            type            surfaceFieldValue;
            libs            ("libfieldFunctionObjects.so");

            writeControl    timeStep;
            writeInterval   ''' + str(int(writeInterval)) +''';

            writeFields     false;	
            operation       areaIntegrate;
            fields          (heatfluxfields);
            regionType      patch;
            name            Outlet;
            region          fluid;
        }
    }

    outletflowrate
    {
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");

        patch           Outlet;
        fields          (phi);

        operation       sum;
        executeControl  timeStep;
        executeInterval ''' + str(int(writeInterval)) +''';
        writeControl    timeStep;
        writeInterval   ''' + str(int(writeInterval)) +''';

        region          fluid;
        writeFields     false;
        regionType      patch;
        name            $patch;
    }
}

// ************************************************************************* //
'''
    f.write(Control_Dict)
    f.close()

def write_fluid_fvSchemes(simulation_name):
    f = open(simulation_name + "/system/fluid/fvSchemes", "w")
    fluid_fvSchemes = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;

    div(phi,U)      bounded Gauss upwind;
    div(phi,K)      bounded Gauss upwind;
    div(phi,h)      bounded Gauss upwind;
    div(phi,k)      bounded Gauss upwind;
    div(phi,K)      bounded Gauss upwind;
    div(phi,R)      bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div(R)          Gauss linear;
    div(Ji,Ii_h)    Gauss linearUpwind grad(U);
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default        Gauss linear uncorrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         uncorrected;
}

wallDist
{
    method meshWave;
}

// ************************************************************************* //
'''
    f.write(fluid_fvSchemes)
    f.close()

def write_fluid_fvSolution(simulation_name):
    f = open(simulation_name + "/system/fluid/fvSolution", "w")
    fluid_fvSolution = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    rho
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0;
    }

    p_rgh
    {
        solver           GAMG;
        tolerance        1e-7;
        relTol           0.01;

        smoother         GaussSeidel;

    }

    "(U|h|e|k|G|Ii)"
    {
        solver           PBiCGStab;
        preconditioner   DILU;
        tolerance        1e-7;
        relTol           0.1;
    }
    omega
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
        nSweeps         1;
    }

    G
    {
        $p_rgh;
        tolerance       1e-05;
        relTol          0.1;
    }
}

PIMPLE
{
    nNonOrthogonalCorrectors 0;
}

relaxationFactors
{
    fields
    {
        rho             1.0;
        p_rgh           0.7;
    }
    equations
    {
        U               0.3;
        h		0.90;
        "(k|epsilon|omega)" 0.7;
        G               0.7;
        "ILambda.*"     0.7;
        qr              0.7;
    }
}

// ************************************************************************* //
'''
    f.write(fluid_fvSolution)
    f.close()

def write_solid_fvSchemes(simulation_name):
    f = open(simulation_name + "/system/solid/fvSchemes", "w")
    solid_fvSchemes = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default     steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
}

laplacianSchemes
{
    default             none;
    laplacian(alphae,e)  Gauss linear uncorrected;
    laplacian(thermo:kappa,T)  Gauss linear uncorrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         uncorrected;
}


// ************************************************************************* //
'''
    f.write(solid_fvSchemes)
    f.close()

def write_solid_fvOptions(simulation_name):
    f = open(simulation_name + "/system/solid/fvOptions", "w")
    solid_fvOptions = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvOptions;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
'''
    f.write(solid_fvOptions)
    f.close()

def write_solid_fvSolution(simulation_name):
    f = open(simulation_name + "/system/solid/fvSolution", "w")
    solid_fvSolution = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    e
    {
        solver           PCG;
        preconditioner   DIC;
        tolerance        1e-06;
        relTol           0.1;
    }
}

PIMPLE
{
    nNonOrthogonalCorrectors 0;
}

relaxationFactors
{
    equations
    {
        e               1.0;
    }
}

// ************************************************************************* //
'''
    f.write(solid_fvSolution)
    f.close()

def write_solid_temperature(simulation_name, Temperature_Inlet, Heat_Flux):
    f = open(simulation_name + "/0/solid/T", "w")
    solid_temperature = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       volScalarField;
    location    "0/solid";
    object      T;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 0 0 1 0 0 0 ];

internalField   uniform ''' + str(Temperature_Inlet) + ''';

boundaryField
{
    Adiabatic_Walls
    {
        type            zeroGradient;
    }
    Heated_Face
    {
        type            externalWallHeatFluxTemperature; //Pour heat source ou Heat Flux
        mode            power;
        Q               ''' + str(Heat_Flux) + ''';
        kappaMethod     solidThermo;
        value           uniform ''' + str(Temperature_Inlet) + ''';
    }
    solid_to_fluid
    {
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        kappaMethod     solidThermo;
        value           uniform ''' + str(Temperature_Inlet) + ''';
    }
}

// ************************************************************************* //
'''
    f.write(solid_temperature)
    f.close()

def write_fluid_velocity(simulation_name, mass_flow_rate, rho):
    f = open(simulation_name + "/0/fluid/U", "w")
    fluid_velocity = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       volVectorField;
    location    "0/fluid";
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 1 -1 0 0 0 0 ];

internalField   uniform ( 0 0 0 );

boundaryField
{
    Inlet
    {
        type            flowRateInletVelocity;
        massFlowRate    ''' + str(mass_flow_rate) + ''';
        rho             rho;
        rhoInlet	''' + str(rho) + ''';
    }
    Outlet
    {
        type            zeroGradient;
    }
    fluid_to_solid
    {
        type            noSlip;
    }
}

// ************************************************************************* //
'''
    f.write(fluid_velocity)
    f.close()

def write_fluid_pressure(simulation_name):
    f = open(simulation_name + "/0/fluid/p", "w")
    fluid_pressure = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       volScalarField;
    location    "0/fluid";
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 1 -1 -2 0 0 0 0 ];

internalField   uniform 0;

boundaryField
{
    Inlet
    {
        type            zeroGradient;
    }
    Outlet
    {
        type            fixedValue;
        value           uniform 0;
    }
    fluid_to_solid
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //
'''
    f.write(fluid_pressure)
    f.close()

def write_fluid_pressure_rgh(simulation_name):
    f = open(simulation_name + "/0/fluid/p_rgh", "w")
    fluid_pressure_rgh = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       volScalarField;
    location    "0/fluid";
    object      p_rgh;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 1 -1 -2 0 0 0 0 ];

internalField   uniform 0;

boundaryField
{
    Inlet
    {
        type            zeroGradient;
    }
    Outlet
    {
        type            fixedValue;
        value           uniform 0;
    }
    fluid_to_solid
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //
'''
    f.write(fluid_pressure_rgh)
    f.close()

def write_fluid_temperature(simulation_name, Temperature_Inlet):
    f = open(simulation_name + "/0/fluid/T", "w")
    fluid_temperature = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       volScalarField;
    location    "0/fluid";
    object      T;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 0 0 1 0 0 0 ];

internalField   uniform ''' + str(Temperature_Inlet) + ''';

boundaryField
{
    Inlet
    {
        type            fixedValue;
        value           uniform ''' + str(Temperature_Inlet) + ''';
    }
    Outlet
    {
        type            zeroGradient;
    }
    fluid_to_solid
    {
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        value           uniform ''' + str(Temperature_Inlet) + ''';
        Tnbr            T;
        kappaMethod     fluidThermo;
    }
}

// ************************************************************************* //
'''
    f.write(fluid_temperature)
    f.close()

def write_alphat(simulation_name):
    f = open(simulation_name + "/0/fluid/alphat", "w")
    alphat = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       volScalarField;
    location    "0/fluid";
    object      alphat;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 1 -1 -1 0 0 0 0 ];

internalField   uniform 0;

boundaryField
{
    Inlet
    {
        type            calculated;
        value           uniform 0;
    }
    Outlet
    {
        type            calculated;
        value           uniform 0;
    }
    fluid_to_solid
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
}

// ************************************************************************* //
'''
    f.write(alphat)
    f.close()
def write_turbulent_intensity(simulation_name, k):
    f = open(simulation_name + "/0/fluid/k", "w")
    turbulent_intensity = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       volScalarField;
    location    "0/fluid";
    object      k;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 2 -2 0 0 0 0 ];

internalField   uniform ''' + str(k) +''';

boundaryField
{
    Inlet
    {
        type            fixedValue;
        value           $internalField;
    }
    Outlet
    {
        type            zeroGradient;
    }
    fluid_to_solid
    {
        type            kqRWallFunction;
        value           $internalField;
    }
}

// ************************************************************************* //
'''
    f.write(turbulent_intensity)
    f.close()

def write_turbulent_viscosity(simulation_name, nut):
    f = open(simulation_name + "/0/fluid/nut", "w")
    turbulent_viscosity = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       volScalarField;
    location    "0/fluid";
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 2 -1 0 0 0 0 ];

internalField   uniform ''' + str(nut) + ''';

boundaryField
{
    Inlet
    {
        type            calculated;
        value           uniform 0;
    }
    Outlet
    {
        type            calculated;
        value           uniform 0;
    }
    fluid_to_solid
    {
        type            nutkWallFunction;
        value           $internalField;
    }
}

// ************************************************************************* //
'''
    f.write(turbulent_viscosity)
    f.close()

def write_turbulent_dissipation_rate(simulation_name, omega, mixing_length):
    f = open(simulation_name + "/0/fluid/omega", "w")
    turbulent_dissipation_rate = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       volScalarField;
    location    "0/fluid";
    object      omega;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

internalField   uniform ''' + str(omega) + ''';

boundaryField
{
    Inlet
    {
        type            turbulentMixingLengthFrequencyInlet;
        mixingLength    ''' + str(mixing_length) + ''';
        value           $internalField; 
    }
    Outlet
    {
        type            zeroGradient;
    }
    fluid_to_solid 
    {
        type            omegaWallFunction;
	value 		$internalField ;
    }
}

// ************************************************************************* //
'''
    f.write(turbulent_dissipation_rate)
    f.close()

def write_density(simulation_name, rho):
    f = open(simulation_name + "/0/fluid/rho", "w")
    density = header + '''
FoamFile
{
    version     10;
    format      ascii;
    class       volScalarField;
    location    "0/fluid";
    object      rho;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -3 0 0 0 0 0];

internalField   uniform ''' + str(rho) + ''';

boundaryField
{
    Inlet
    {
        type            calculated;
        value           uniform ''' + str(rho) + ''';
    }
    fluid_to_solid
    {
        type            calculated;
        value           uniform ''' + str(rho) + ''';
    }
    Outlet
    {
        type            calculated;
        value           uniform ''' + str(rho) + ''';
    }
}


// ************************************************************************* //
'''
    f.write(density)
    f.close()