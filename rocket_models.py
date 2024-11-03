import numpy as np

gamma = 1.4 #TODO: Verify correct adiabatic index for air over range

class OrbitingSatellite:
    def __init__(self):
        self.mass_satellite = 100           # kg
        self.satellite_length = 0.25        # m
        self.satellite_width = 0.25         # m
        self.satellite_height = 0.25        # m
        self.cross_section_area = 0.0625    # m^2

    def getMass(self):
        return self.mass_satellite
    
    def getDragCoefficient(self, air_speed, density, pressure):
        return 1.0
    
    def getGravityCenter(self):
        return 0.0
    
    def getMOIs(self):
        Ixx = 1/12*self.mass_satellite*(self.satellite_length**2 + self.satellite_height**2 )
        Iyy = 1/12*self.mass_satellite*(self.satellite_height**2 + self.satellite_width**2 )
        Izz = 1/12*self.mass_satellite*(self.satellite_length**2 + self.satellite_width**2 )
        return Ixx, Iyy, Izz
    
    def getThrustVector(self, twist, pressure):
        return 0.0, [0.0,0.0,0.0], 0.0
    
    def getCenterofPressure(self):
        return 0.0
    
class FalconIX:
    def __init__(self):
        #TODO: params
        self.gimbalAngleMax = 5.0*np.pi/180                                 # rad
        self.mass_satellite = 100.0                                         # kg
        self.stage_one_fuselage_mass = 25600.0-470.0*9                      # kg - https://www.researchgate.net/publication/363319640_AAS_22-821_AMBIGUITY_REMEDIATION_IN_SPACE_LAUNCH_VEHICLES_WITH_PARAMETER_UNCERTAINTIES_A_COMPARISON_BETWEEN_SPECIAL_EUCLIDEAN_GROUP_AND_DUAL_QUATERNIONS
        stage_one_liquid_mass = 92670.0                                     # kg
        self.stage_one_fuel_mass = stage_one_liquid_mass*2.6/3.6            # kg
        self.stage_one_oxidizer_mass = stage_one_liquid_mass*1.0/3.6        # kg
        self.stage_one_motor_mass = 470.0*9                                 # kg - Merlin 1D - https://en.wikipedia.org/wiki/SpaceX_Merlin
        self.stage_two_fuselage_mass = 2900.0 - 470.0                       # kg
        stage_two_liquid_mass = 395700+19600                                # kg
        self.stage_two_fuel_mass = stage_two_liquid_mass*2.6/3.6            # kg
        self.stage_two_oxidizer_mass = stage_two_liquid_mass*1.0/3.6        # kg
        self.stage_two_motor_mass = 470.0                                   # kg

        self.stage_one_exit_velocity = 3000.0                               # m/s   9 Merlin 1D Motors
        self.stage_one_exit_pressure = 70927.5                              # N/m^2
        self.stage_one_exit_area = 6.1316                                   # m^2
        self.stage_one_fuel_consumption_rate = 2100.0                       # kg/s

        self.stage_two_exit_velocity = 3000.0                               # m/s - 1 Merlin 1D Motor
        self.stage_two_exit_pressure = 70927.5                              # N/m^2
        self.stage_two_exit_area = 0.68128896                               # m^2
        self.stage_two_fuel_consumption_rate = 2100.0/9                     # kg/s   

        # Dimensional References: https://www.researchgate.net/figure/Falcon-9-launch-vehicle-dimensions-and-mass-breakdown_fig1_363319640
        reference_dimension = 58.0                                          # m
        merlin_length = 2.92                                                # m https://en.wikipedia.org/wiki/SpaceX_Merlin
        self.stage_one_motor_distance = reference_dimension - 0.0           # m
        self.stage_one_oxidizer_max_distance = reference_dimension - 0.384  # m
        self.stage_one_fuel_max_distance = reference_dimension - 17.8       # m
        self.stage_two_motor_distance = reference_dimension - 47.5 + merlin_length  # m
        self.stage_two_oxidizer_max_distance = reference_dimension - 47.5   # m
        self.stage_two_fuel_max_distance = reference_dimension - 51.0       # m

        self.cross_section_area = np.pi*5.2**2                              #m^2
        self.rocket_outer_diameter = 3.65                                   # m
        stage_one_mass = self.stage_one_fuselage_mass+stage_one_liquid_mass+self.stage_one_motor_mass
        stage_two_mass = self.stage_two_fuselage_mass+stage_two_liquid_mass+self.stage_two_motor_mass
        self.total_mass = self.mass_satellite+stage_one_mass+stage_two_mass

        # Assumption - fuselage mass is uniform cylinder of outer shell radius as listed, 
        aluminum_lithium_alloy_density = 2700.0                             # kg/m^3 - https://www.wixsteel.com/products/aluminum-alloy/2000-series-aluminum-alloy/2198
        self.inner_radius = np.sqrt(self.rocket_outer_diameter**2 - (self.stage_one_fuselage_mass+self.stage_two_fuselage_mass)/(aluminum_lithium_alloy_density*reference_dimension*np.pi))
        self.stage_flag = 1

    def getDragCoefficient(self,air_speed,density,pressure):
        #ideal gas - source https://en.wikipedia.org/wiki/Speed_of_sound
        
        if density > 0:
            speed_of_sound = np.sqrt(gamma*pressure/density) 
        else:
            speed_of_sound = 0.0001 #Arbitrarily Low Speed of Sound in Space

        M = air_speed/speed_of_sound
        # Drag Coefficient Model Source: http://www.braeunig.us/space/aerodyn_wip.htm
        if M <= 0.6:
            Cd = 0.2083333*M**2 - 0.25*M + 0.46
        elif 0.6 < M <= 0.8:
            Cd = 1.25*M**3 -2.125*M**2 + 1.2*M + 0.16
        elif 0.8 < M <= 0.95:
            Cd = 10.37037*M**3 - 22.88889*M**2 + 16.9111*M - 3.78963
        elif 0.95 < M <= 1.05:
            Cd = -30*M**3 + 88.5*M**2 - 85.425*M + 27.51375
        elif 1.05 < M <= 1.15:
            Cd = -20*M**3 + 60*M**2 - 58.65*M + 19.245
        elif 1.15 < M <= 1.3:
            Cd = 11.85185*M**3 - 44.88889*M**2 + 56.22222*M - 22.58519
        elif 1.3 < M <= 2:
            Cd = -0.04373178*M**3 + 0.3236152*M**2 - 1.019679*M + 1.554752
        elif 2 < M <= 3.25:
            Cd = 0.01024*M**3 - 0.00864*M**2 - 0.33832*M + 1.08928
        elif 3.25 < M <= 4.5:
            Cd = -0.01408*M**3 + 0.19168*M**2 - 0.86976*M + 1.53544
        else:
            Cd = 0.22

        return Cd
    
    def getGravityCenter(self):
        # Considered relative to payload Cg
        # Cg is on centerline of rocket

        if self.stage_flag == 1:
            Cg = 0.0
        if self.stage_flag == 2:
            Cg = 0.0
        else:
            Cg = 0.0

        return Cg
    
    def getMOIs(self):
        Ixx = 2945907.90    #stolen from Saturn V - only inital conditions
        Iyy = 892137606
        Izz = 892137606

        return Ixx, Iyy, Izz

    def getThrustVector(self, twist, air_pressure):
        #TODO: Thrust Vectoring Functions
        # Thrust Calcs - https://space.stackexchange.com/questions/46521/falcon-9-merlin-1d-thrust-calculated-through-every-moment-of-flight
        if self.stage_flag == 1:
            F_thrust = self.stage_one_fuel_consumption_rate*self.stage_one_exit_velocity+(self.stage_one_exit_pressure-air_pressure)*self.stage_one_exit_area
            self.stage_one_fuel_mass -= self.stage_one_fuel_consumption_rate*2.66/3.66
            self.stage_one_oxidizer_mass -= self.stage_one_fuel_consumption_rate*1.00/3.66

            motor_distance = self.stage_one_motor_distance
            if self.stage_one_fuel_mass <= 0: # There is a slight amount of error here equal to the negative magnitude of fuel used TODO: incorporate stage separation properly
                self.stage_flag = 2
        

        elif self.stage_flag == 2:
            F_thrust = self.stage_two_fuel_consumption_rate*self.stage_two_exit_velocity+(self.stage_two_exit_pressure-air_pressure)*self.stage_two_exit_area
            self.stage_two_fuel_mass -= self.stage_two_fuel_consumption_rate*2.66/3.66
            self.stage_two_oxidizer_mass -= self.stage_two_fuel_consumption_rate*1.00/3.66
            motor_distance = self.stage_two_motor_distance
            if self.stage_two_fuel_mass < 0:
                self.stage_flag = 3
        
        else:
            F_thrust = 0.0
            motor_distance = 0.0

        # TODO: Develop Gimbal Control Method and Control Law (also Trajectories)
        return F_thrust, [1.0, 0.0, 0.0], motor_distance
    
    def getCenterofPressure(self):
        if self.stage_flag == 1:
            return 0.0
        elif self.stage_flag == 2:
            return 0.0
        else:
            return 0.0
        
    def getMass(self):
        return self.total_mass