import numpy as np

gamma = 1.4 #TODO: Verify correct adiabatic index for air over range

class OrbitingSatellite:
    def __init__(self):
        self.mass_satellite = 100           #kg

    def getMass(self):
        return(self.mass_satellite)
    
    def getDragCoefficient(self,air_speed, density, pressure):
        return(0.0)
    
    def getGravityCenter(self):
        return(0.0)
    
class LiquidRocket:
    def __init_(self):
        #TODO: params
        self.gimbalAngleMax = 4.5           # deg
        self.mass_satellite = 100           # kg
        self.mass_booster = 0.0             # kg
        self.mass_booster_fuel = 0.0        # kg
        self.mass_second_stage = 0.0        # kg
        self.mass_second_stage_fuel = 0.0   # kg
        # Center of lift and dimensions source:

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
        return(0.0)

    def getThrustVector(self):
        #TODO: Thrust Vectoring Functions
        return None