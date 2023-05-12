''' AIRCRAFT MASSES & LOCATIONS OUTPUTS --> ALBERT  ---> Last Update - 08/03/23 '''

''' Overall system masses '''
m_f = 636.7765587973576
trip_fuel = 400 #419.97377851582326
reserve_fuel = m_f - trip_fuel #242.45814995926958
m_oem = 12634.92498884976
m_mto = 18777.90521758506
m_zf = 18141.1286587877
eta_eng = 0.6 * 0.97 * 0.995 * 0.95
m_pldes = 5524.75056
W_S_design = 2820.64
W_P_design = 0.07816348375747666
beta_s_land_fc = m_zf / m_mto

''' Subsystem masses '''

''' Sub system locations '''

''' Overall CG locations '''
LEMAC = 12.087743576790354  # LEMAC position                    [m]
xcg_gear = 13.4             # Location landing gear             [m]
xcg_aft_potato = 13.065460618085575       # Aft cg location (potato plot)     [m]
xcg_front_potato = 12.21880418414821     # Front cg location (potato plot)   [m]
cg_tank = 1.3
