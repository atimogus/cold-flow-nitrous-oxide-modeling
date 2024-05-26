from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

temp_degrees = 0.5
oxidizer_mass_l_kg =19
time_step = 0.03
fluid , fuel = 'N2O' , 'Acetone'
o_f_ratio = 3
P_atmosphere = 1e5
liquid , gas = 0 , 1
gamma = 1.3
R = PropsSI('GAS_CONSTANT', 'NitrousOxide') / PropsSI('M', 'NitrousOxide')  # Specific gas constant in J/(kg*K)
temp = 273.15 + (temp_degrees) #pressure of a tank is dependant on a tank temperature

density_l = PropsSI('D','T',temp,'Q',liquid,fluid)
V_liquid = oxidizer_mass_l_kg/density_l  # additional 10 percent for gas
V_gas =0.01*V_liquid
V_tank = V_liquid + V_gas
mass_v = V_gas * PropsSI('D','T',temp,'Q',gas,fluid)
mass_l = oxidizer_mass_l_kg

#oxidizer injector specs
C_d = 0.44
d_mm = 7.2
A = (d_mm/1000)**2 * np.pi/4
pressure_l = PropsSI('P','T',temp,'Q',liquid,fluid)

#fuel specs & fuel injector specs
density_f = PropsSI('D', 'P', pressure_l, 'T', temp, fuel)
mass_f = oxidizer_mass_l_kg / (o_f_ratio * 1.5)
mass_f_initial = mass_f
V_fuel = mass_f / density_f

pressure_frictionLoss  = 2e5 #needs to be determined experimentaly
ml_dot_initial = C_d * A * np.sqrt(2 * density_l * (pressure_l - P_atmosphere))
A_fuel = ml_dot_initial / (o_f_ratio *(C_d * np.sqrt(2 * density_f * (pressure_l - pressure_frictionLoss - P_atmosphere)))) #guessing fuel size orfice based on mass flow rate of oxidizer

liquid_pressure = []
liquid_massflow = []
vapor_pressure = []
vapor_mass_flow = []
mass_liquid = []
mass_vapour = []
O_F_shift = []
fuel_massflow = []
fuel_pressure = []
dp_dt = []

while mass_l > oxidizer_mass_l_kg * 0.01:
    mass_liquid , mass_vapour =  np.append(mass_liquid, mass_l) , np.append(mass_vapour, mass_v)
    
    pressure_old = pressure_l
    density_v, pressure_v = PropsSI('D','T',temp,'Q',gas,fluid), PropsSI('P','T',temp,'Q',gas,fluid) # promjena gustoce plina
    density_l, pressure_l = PropsSI('D','T',temp,'Q',liquid,fluid), PropsSI('P','T',temp,'Q',liquid,fluid)# promjena gustoce tekucine, pressure vapour = pressure liquid

    if  mass_l < oxidizer_mass_l_kg *0.9  and pressure_old - pressure_l < np.mean(dp_dt)*0.5:
        temp = PropsSI('T','P',pressure_old - np.mean(dp_dt),'Q',gas,fluid)
        density_v, pressure_v = PropsSI('D','T',temp,'Q',gas,fluid), PropsSI('P','T',temp,'Q',gas,fluid) # promjena gustoce plina
        density_l, pressure_l = PropsSI('D','T',temp,'Q',liquid,fluid), PropsSI('P','T',temp,'Q',liquid,fluid)# promjena gustoce tekucine, pressure vapour = pressure liquid
    dp_dt = np.append(dp_dt, pressure_old - pressure_l)

    liquid_pressure = np.append(liquid_pressure, pressure_v)
    #dizajn injektora za jednu rupu , treba promjeniti za odgovorarajuce dizajne
    ml_dot = C_d * A * np.sqrt(2 * density_l * (pressure_l - P_atmosphere)) * time_step
    liquid_massflow = np.append(liquid_massflow, ml_dot)
    mass_l -= ml_dot 
    mass_l_old = mass_l
    
    if mass_f > mass_f_initial * 0.01:
        pressure_f = pressure_l - pressure_frictionLoss
        density_f = PropsSI('D', 'P', pressure_l, 'T', temp , fuel)
        mf_dot = C_d * A_fuel * np.sqrt(2 * density_f * (pressure_f - P_atmosphere)) * time_step
        mass_f -= mf_dot
        
        V_tank +=  mf_dot / density_f
        mass_l = (V_tank-(mass_l_old+mass_v)/density_v)/(1/density_l-1/density_v)
        mass_vapourised = mass_l_old - mass_l
        mass_vapourised_old = mass_vapourised
        mass_v += mass_vapourised

        fuel_pressure = np.append(fuel_pressure, pressure_f)
        O_F_shift = np.append(O_F_shift, ml_dot/mf_dot)
        fuel_massflow = np.append(fuel_massflow, mf_dot)
        
    else:
        mf_dot = 0
        fuel_massflow, fuel_pressure = np.append(fuel_massflow, mf_dot), np.append(fuel_pressure, mf_dot)
        O_F_shift = np.append(O_F_shift, mf_dot)

        mass_l = (V_tank-(mass_l_old+mass_v)/density_v)/(1/density_l - 1/density_v)        
        mass_vapourised = (mass_l_old - mass_l)
        mass_v += mass_vapourised
    
    Heat_of_vapourisation = PropsSI('H','T',temp,'Q',gas,fluid) - PropsSI('H','T',temp,'Q',liquid,fluid)
    heat_removed_deltaQ = mass_vapourised * Heat_of_vapourisation # entalpija tekucine - entalpija plina = PropsSI('H','T',300,'Q',0,'N2O') - PropsSI('H','T',300,'Q',1,'N2O')
    deltaT = -heat_removed_deltaQ/(mass_l * PropsSI('C','T',temp,'Q', liquid ,fluid))
    temp += deltaT #promjena temperature tekucine usljed hladjenja tekucine zbog isparavanja       

# Initial settings for vapor loop
mass_v_initial = mass_v
# Modeling for vapor pressure as ideal gas
while pressure_v > P_atmosphere * 1 and mass_v > mass_v_initial * 0.05: 
    fuel_massflow, fuel_pressure = np.append(fuel_massflow, mf_dot), np.append(fuel_pressure, mf_dot)
    mass_liquid , mass_vapour =  np.append(mass_liquid, mass_l) , np.append(mass_vapour, mass_v)
    O_F_shift = np.append(O_F_shift, mf_dot)

    mv_dot = C_d * A * np.sqrt(2 * density_v * (pressure_v - P_atmosphere)) * time_step
    mass_old = mass_v
    mass_v -= mv_dot

    # Update temperature and pressure using polytropic relations
    T_2 = temp * (mass_v / mass_old) ** (gamma - 1)
    P_2 = pressure_v * (T_2 / temp) ** (gamma / (gamma - 1))
    temp = T_2
    pressure_v = P_2

    vapor_pressure.append(P_2)
    vapor_mass_flow.append(mv_dot)

tank_pressure = np.concatenate((liquid_pressure , vapor_pressure))
tank_massflow = np.concatenate((liquid_massflow, vapor_mass_flow)) / time_step
fuel_massflow = fuel_massflow / time_step
sequence = np.arange(0, len(tank_massflow)*time_step,time_step)  

# # Write to CSV file
# with open('tank_data.csv', mode='w', newline='') as file:
#     writer = csv.writer(file)
#     writer.writerow(['Time (s)', 'Tank Pressure (Pa)', 'Mass Flow Rate (kg/s)'])
#     for i in range(len(sequence)):
#         writer.writerow([sequence[i], tank_pressure[i], tank_massflow[i]])

dataset = pd.read_csv('Dataset.csv')

# Plotting results
fig = 6
plt.figure(figsize=(fig * 10, fig))

# Plot for vapor tank pressure
plt.subplot(1, 3, 1)
plt.plot(sequence, tank_pressure/1e5, color='red')
plt.plot(dataset['timeShift'], dataset['pressureShift']-3, label='Dataset Pressure', color='blue')
plt.plot(sequence, fuel_pressure/1e5, color='green')
plt.xlabel('Time (s)')
plt.ylabel('tank pressure (Pa)')
plt.legend(['ox cold flow model', 'experimantal cold flow data', 'fuel cold flow model'])
plt.grid(True)

# Plot for vapor mass flow
plt.subplot(1, 3, 2)
plt.plot(sequence, tank_massflow, color='red')
plt.plot(sequence, fuel_massflow, color='green')
plt.xlabel('Time (s)')
plt.ylabel('massflow (kg/s)')
plt.legend(['N2O mass flow rate', 'Fuel mass flow rate'])
plt.grid(True)

# Plot for vapor tank pressure
plt.subplot(1, 3, 3)
plt.plot(sequence, mass_liquid, color='red')
plt.plot(sequence, mass_vapour, color='green')
plt.xlabel('Time (s)')
plt.ylabel('mass in tank (kg)')
plt.legend(['mass_liquid', 'mass_vapour'])
plt.grid(True)

plt.tight_layout()
plt.show()

print(f"mean ox mass flow = {np.mean(liquid_massflow) / time_step:.2f}, mean fuel mass flow = {np.mean(fuel_massflow[fuel_massflow != 0]):.2f}")
