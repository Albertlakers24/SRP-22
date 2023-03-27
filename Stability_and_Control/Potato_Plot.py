from Constants.MissionInputs import m_pax_baggage,vol_pax_baggage,PAX,m_pax
from matplotlib import pyplot as plt
import numpy as np
from Constants.Masses_Locations import m_oem, m_f
from Constants.AircraftGeometry import c_mac_w

cg_start = 12.33#475                                    #OEW center of gravity from nose (m)
cargo_front = 123                                       #Front cargo position from the nose (m)
cargo_aft = 13                                          #Aft cargo position from the nose (m)
l_cabin = 13.0                                          #Cabin length (m)
overhead_surface = 0.1508                               #Surface area per overhead bin (m^2)
overhead_volume = overhead_surface * 2 * l_cabin        #Overhead volume available
req_volume = vol_pax_baggage                            #Cargo volume required
front_volume = max(0, req_volume - overhead_volume)     #Cargo volume in front
print(front_volume, "VOLUME FRONT", overhead_volume, "VOLUME OVERHEAD")
front_surface_area = 1.37                               #Cargo surface area per front compartment (m^2)
len_front_storage = front_volume/(2*front_surface_area) #Length of front compartment (m)
mass_front = front_volume/vol_pax_baggage*m_pax_baggage #Mass front storage compartment (kg)
mass_overhead = m_pax_baggage - mass_front              #Mass overhead storage (kg)
x_front = 3.3                                           #Start of circular part (m)
x_cg_front = 3.3 + len_front_storage / 2                #CG location of front storage (m)
#POINT LOADS APPLY AT EVERY SEAT, SO 13/12, / 2 to obtain its cg
x_cg_bags = l_cabin / (PAX / 4) / 2
overhead_incr = mass_overhead / 12
x_cabin_start = x_front + len_front_storage
x_cg_over_front = (cg_start + x_cabin_start) / 2
longer_len_cg = cg_start #+ len_front_storage
def cg_shift(old_cg, old_mass, added_cg, added_mass):
    new_cg = (old_cg * old_mass + added_cg * added_mass) / (old_mass + added_mass)
    total_mass = old_mass + added_mass
    return new_cg, total_mass
front_fraction = (longer_len_cg - x_cabin_start) / 13
front_mass = front_fraction * mass_overhead
front_cargo_cg, front_cargo_m = cg_shift(x_cg_front, mass_front, x_cg_over_front, front_mass)
AC_cargo_f_cg, AC_cargo_f_m = cg_shift(longer_len_cg, m_oem, front_cargo_cg, front_cargo_m)
mass_aft = (1 - front_fraction) * mass_overhead
cg_aft = (x_cabin_start * 2 + 13) / 2
AC_cargo_aft_cg, AC_cargo_aft_m = cg_shift(longer_len_cg, m_oem, cg_aft, mass_aft)
AC_cargo_cg_1, AC_cargo_m_1 = cg_shift(AC_cargo_f_cg, AC_cargo_f_m, cg_aft, mass_aft)
AC_cargo_cg_2, AC_cargo_m_2 = cg_shift(AC_cargo_aft_cg, AC_cargo_aft_m, front_cargo_cg, front_cargo_m)
x_vals = [longer_len_cg, AC_cargo_f_cg, AC_cargo_aft_cg, AC_cargo_cg_1, AC_cargo_cg_2]
y_vals = [m_oem, AC_cargo_f_m, AC_cargo_aft_m, AC_cargo_m_1, AC_cargo_m_2]
loading_pax = 2 * m_pax / PAX
seat_incr = l_cabin / PAX * 4

#LOAD FUEL ONLY
def fuel_only():
    tank_mass = m_f
    cg_tank_total = x_cabin_start + l_cabin + cg_tank
    cg_fuel, mass_fuel = cg_shift(cg_start, 12770, cg_tank_total, tank_mass)
    return cg_fuel, mass_fuel

def seating(where, start_cg, start_mass):
    cgs = start_cg#[x_vals[3]]
    masses = start_mass#[y_vals[3]]
    if where == "front":
        for i in range(0, 12):
            if i == 0:
                x_cg = x_cabin_start + 0.5 * seat_incr * (i + 0.5)
                x_cg_new, mass_new = cg_shift(cgs[i], masses[i], x_cg, loading_pax)
            else:
                x_cg = x_cabin_start + seat_incr * (i + 0.5)
                x_cg_new, mass_new = cg_shift(cgs[i], masses[i], x_cg, loading_pax)
            cgs.append(x_cg_new)
            masses.append(mass_new)
    elif where == "back":
        for i in range(0, 12):
            if i == 0:
                x_cg = x_cabin_start + (l_cabin - seat_incr * (i + 0.5))
                x_cg_new, mass_new = cg_shift(cgs[i], masses[i], x_cg, loading_pax)
            else:
                x_cg = x_cabin_start + (l_cabin - seat_incr * (i + 0.5))
                x_cg_new, mass_new = cg_shift(cgs[i], masses[i], x_cg, loading_pax)
            cgs.append(x_cg_new)
            masses.append(mass_new)
    return cgs, masses

cgs_1, masses_1 = seating("front", [x_vals[3]], [y_vals[3]])
cgs_2, masses_2 = seating("back", [x_vals[3]], [y_vals[3]])
cgs_3, masses_3 = seating("front", [cgs_1[12]], [masses_1[12]])
cgs_4, masses_4 = seating("back", [cgs_1[12]], [masses_1[12]])

#LOAD FUEL
tank_mass = m_f #*2.4
cg_tank_total = x_cabin_start + l_cabin + cg_tank
cg_end, mass_end = cg_shift(cgs_3[12], masses_3[12], cg_tank_total, tank_mass)

LEMAC =11.507
c_bar = c_mac_w
def Location_in_MAC(Location):
    """
    :return: Location (MAC)
    """
    xcg_MAC = (Location-LEMAC)/c_bar
    return xcg_MAC
print("OLD VALS", x_vals)
for i in np.arange(0, len(x_vals)):
    x_vals[i] = Location_in_MAC(x_vals[i])
print("NEW VALS", x_vals)
for i in np.arange(0, len(cgs_1)):
    cgs_1[i] = Location_in_MAC(cgs_1[i])
for i in np.arange(0, len(cgs_2)):
    cgs_2[i] = Location_in_MAC(cgs_2[i])
for i in np.arange(0, len(cgs_3)):
    cgs_3[i] = Location_in_MAC(cgs_3[i])
for i in np.arange(0, len(cgs_4)):
    cgs_4[i] = Location_in_MAC(cgs_4[i])
cg_end = Location_in_MAC(cg_end)
#PROPER FUEL
cg_fuel_only = Location_in_MAC(fuel_only()[0])


plt.figure()
plt.plot([Location_in_MAC(cg_start), cg_fuel_only], [m_oem, fuel_only()[1]], "hotpink", label="Fuel First", marker="s")
plt.plot([0, cg_fuel_only], [fuel_only()[1], fuel_only()[1]], 'grey', linestyle='--')
plt.plot([x_vals[0], x_vals[1]], [y_vals[0], y_vals[1]], "black", label="Cargo")
plt.plot([x_vals[0], x_vals[2]], [y_vals[0], y_vals[2]], "black")
plt.plot([x_vals[1], x_vals[3]], [y_vals[1], y_vals[3]], "black")
plt.plot([x_vals[2], x_vals[4]], [y_vals[2], y_vals[4]], "black")
plt.plot(cgs_1, masses_1, "red", marker = "s", label = "Window Seats")
plt.plot(cgs_2, masses_2, "green", marker = "s")
plt.plot(cgs_3, masses_3, "blue", marker = "s", label = "Aisle Seats")
plt.plot(cgs_4, masses_4, "orange", marker = "s")
plt.plot([cgs_3[12], cg_end], [masses_3[12], mass_end], "olive", marker = "s", label = "Fuel")
plt.plot([min(cgs_3) - 0.1, cg_end], [mass_end, mass_end], color="grey", linestyle="--")
plt.plot([min(cgs_3) - 0.1, cgs_3[12]], [masses_3[12], masses_3[12]], color = 'grey', linestyle="--")
plt.annotate('MTOM', xy=(11.5, mass_end + 100))
plt.annotate('MZFM', xy=(11.5, masses_3[12] + 100))
plt.annotate('m_oem + Max Fuel', xy=(11.5, fuel_only()[1] + 100))
plt.xlabel("xcg [MAC]")
plt.ylabel("Mass (kg)")
plt.ylim(m_oem, 19000)
plt.xlim(min(cgs_3) - 0.05, cg_fuel_only + 0.05)
plt.legend()
plt.show()
print(cgs_1)
print(cgs_2)
print(masses_1)
print(masses_2)
print(len_front_storage)
#PRINT AFT AND FRONT CG
print(cg_fuel_only, cgs_3[12])