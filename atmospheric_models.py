r_earth = 6378137 #m
G = 6.67430*10**-11# Earth Parameters - Point Mass Model
g0 = 9.80665  # m/s^2
M = 0.0289644 # kg/mol
R = 8.3144598 # J/molK
e = 2.718281828459
import numpy as np
import bisect

# Standard Atmosphere Source: https://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
class US_Standard_Atmosphere:
    def __init__(self):
        altitudes_from_sea_level = np.array([0,11000,20000,32000,47000,51000,71000,84852])
        self.altitudes = altitudes_from_sea_level+r_earth
        self.temperatures_K = [288.15,216.65,216.65,228.65,270.65,270.65,214.65,186.946]
        self.temperature_lapse_rates = [0.0065,0.0,-0.001,-0.0028,0.0,0.0028,0.002,0.0] #Kelvin per meter
        self.pressures_Pa = [101325,22632.1,5474.89,868.02,110.91,66.94,3.96,0.0]
        self.pressure_exponents = [5.25588, None, -34.1626, -12.2009, None, 12.2009, 17.0813,0.0]

    def get_params(self, altitude):
        temperature = np.interp(altitude, self.altitudes,self.temperatures_K)
        
        # Finds index using floor method to determine which atmospheric range of the US model to reference for calculations
        ref_index = bisect.bisect_left(self.altitudes,altitude)-1
        ref_altitude = self.altitudes[ref_index]
        ref_temperature = self.temperatures_K[ref_index]
        ref_pressure = self.pressures_Pa[ref_index]
        ref_temperature_lapse_rate = self.temperature_lapse_rates[ref_index]
        ref_pressure_exponent = self.pressure_exponents[ref_index]

        if ref_index is not len(self.altitudes)-1:
            if ref_pressure_exponent is not None:
                pressure = ref_pressure*(1-(ref_temperature_lapse_rate/ref_temperature)*(altitude-ref_altitude))**ref_pressure_exponent
            else:
                pressure = ref_pressure*np.exp(-g0*M*(altitude-ref_altitude)/(R*ref_temperature))
        else:
            pressure = 0.0
            temperature = self.temperatures_K[-1]

        #Density (Ideal) Gas Law Source: https://en.wikipedia.org/wiki/Density_of_air
        density = pressure*M/(R*temperature)

        return pressure, temperature, density