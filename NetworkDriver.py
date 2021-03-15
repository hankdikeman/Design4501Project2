"""
This file implements the reactor network optimization algorithm for the OME3-5
production process using the ReactorNetwork object and component objects.
"""
import numpy as np
from ReactorNetwork import Network
from Components import Component, HeatExchanger, Compressor, Turbine, Reactor, DistillationColumn, Absorber, Adsorber, Feed, Removal, ProductRemoval, Mixer, Splitter, FlashTank, StreamGen
from ConversionFunctions import MethanolReactor, FormaldehydeReactor, OMEReactor
import warnings
# ignore annoying numpy truedivide warnings, implement our own for divergence
warnings.filterwarnings("ignore")


if __name__ == "__main__":
    # print docstring
    print(__doc__)
    # make network object
    net1 = Network()
    # feed to reactor
    i1outlets = {'i1': StreamGen(H2=60, CO2=20)}
    net1.add_component('I1', Feed(i1outlets))
    # mixer for recycle and feed
    mix1inlets = {
        'i1': StreamGen(H2=60, CO2=20),
        's11': StreamGen(H2=60, CO2=20, CO=5)
    }
    mix1outlets = {'m11': StreamGen()}
    net1.add_component('M1', Mixer(mix1inlets, mix1outlets))
    # methanol reactor
    TEMP_MEOH = 673
    PRESS_MEOH = 37500
    r1inlets = {'m11': StreamGen(H2=100, CO2=35, MEOH = 5, H2O=5, CO=5)}
    r1outlets = {'r11': StreamGen()}
    net1.add_component('R1', Reactor(TEMP_MEOH, PRESS_MEOH,
                                     r1inlets, r1outlets, MethanolReactor))
    # splitter for recycle
    s1inlets = {'f11': StreamGen(CO2=20, H2=60, CO=5)}
    s1outlets = {
        's11': StreamGen(),
        's12': StreamGen()
    }
    net1.add_component('S1', Splitter(s1inlets, s1outlets, 0.98, 's11'))
    # purge outlet
    purge1inlets = {'s12': StreamGen(CO2=1, CO=1, H2=1)}
    net1.add_component('Purge1', Removal(purge1inlets))
    # flash tank for reactor outlet
    f1inlets = {'r11': StreamGen(CO2=20, H2=60, MEOH=15, H2O=15, CO=5)}
    f1vaporoutlets = {'f11': StreamGen()}
    f1liquidoutlets = {'f12': StreamGen()}
    f1recovery = {'Key': (0.0005, 4)}
    F1PRESS = 760
    net1.add_component('F1', FlashTank(F1PRESS, f1recovery,
                                       4, f1inlets, f1vaporoutlets, f1liquidoutlets))
    # distillation column for methanol
    DC1_recov = {'HK':(0.01,3),'LK':(0.99,4)}
    inlets = {'f12': StreamGen(MEOH=15, H2O=15)}
    v_out = {'dc11': StreamGen()}
    l_out = {'dc12': StreamGen()}
    DC1_TEMP = 350
    DC1_PRESS = 760 * 2
    net1.add_component('DC1', DistillationColumn(DC1_TEMP, DC1_PRESS, DC1_recov, inlets, v_out, l_out))
    # water outlet
    waterinlet = {'dc12': StreamGen(H2O=15)}
    net1.add_component('WaterR1', Removal(waterinlet))
    # splitter for pure methanol
    s2inlets = {'dc11': StreamGen(MEOH=15)}
    s2outlets = {
        's21': StreamGen(),
        's22': StreamGen()
    }
    net1.add_component('S2', Splitter(s2inlets, s2outlets, 0.55, 's21'))
    # air inlet
    i2outlets = {'i2': StreamGen(O2=10, N2=40)}
    net1.add_component('I2', Feed(i2outlets))
    # mixer for oxygen and methanol into formaldehyde
    mix2inlets = {
        's21': StreamGen(MEOH=10),
        'i2': StreamGen(O2=60, N2=20)
    }
    mix2outlets = {'m21': StreamGen()}
    net1.add_component('M2', Mixer(mix2inlets, mix2outlets))
    # formaldehyde reactor
    r2inlets = {'m21': StreamGen(MEOH=10, O2=10, N2=40)}
    r2outlets = {'r21': StreamGen()}
    TEMP_FA = 700
    PRESS_FA = 760
    net1.add_component('R2', Reactor(TEMP_FA, PRESS_FA,
                                     r2inlets, r2outlets, FormaldehydeReactor))
    # water inlet
    i3outlets = {'i3': StreamGen(H2O=2)}
    net1.add_component('I3', Feed(i3outlets))
    # formaldehyde absorber
    vinabs = {'r21': StreamGen(FA=10, O2=5, N2=40, H2O=10)}
    linabs = {'i3': StreamGen(H2O=5)}
    voutabs = {'a11': StreamGen(O2=5, N2=40)}
    loutabs = {'a12': StreamGen()}
    ABS_WATERTEMP = 310
    ABS_PRESS = 760
    abs_key_recov = 0.99
    abs_key_index = 5
    abs_solvent_index = 3
    net1.add_component('A1', Absorber(ABS_WATERTEMP, ABS_PRESS, abs_key_recov, abs_key_index, abs_solvent_index, vinabs, voutabs, linabs, loutabs))
    # air removal
    airinlet  = {'a11': StreamGen()}
    net1.add_component('AirR1', Removal(airinlet))
    # OME + methanol mixer
    mix3inlets = {
        'a12': StreamGen(FA=10, H2O=15),
        'dc21': StreamGen(MEOH=10, FA=18, H2O=0.1, OME1=0.1, OME2=0.1, OME3=0.1, OME4=0.1, OME5=0.1, OME6=0.1)
    }
    mix3outlets = {'m31': StreamGen()}
    net1.add_component('M3', Mixer(mix3inlets, mix3outlets))
    # recycle + OME + methanol mixer
    mix4inlets = {
        'ad11': StreamGen(MEOH=10, FA=18, H2O=0.1, OME1=0.1, OME2=0.1, OME3=0.1, OME4=0.1, OME5=0.1, OME6=0.1),
        's22': StreamGen(MEOH=5)
    }
    mix4outlets = {'m41': StreamGen()}
    net1.add_component('M4', Mixer(mix4inlets, mix4outlets))
    # OME reactor
    r3inlets = {'m41': StreamGen(MEOH=10, FA=18, H2O=0.1, OME1=0.1, OME2=0.1, OME3=0.1, OME4=0.1, OME5=0.1, OME6=0.1)}
    r3outlets = {'r31': StreamGen()}
    TEMP_OME = 350
    PRESS_OME = 760
    net1.add_component('R3', Reactor(TEMP_OME, PRESS_OME,
                                     r3inlets, r3outlets, OMEReactor))
    # OME column
    DC2_recov = {'HK':(0.01,10),'LK':(0.99,9)}
    inletsdc2 = {'r31': StreamGen(MEOH=10, FA=18, H2O=0.1, OME1=0.1, OME2=0.1, OME3=0.1, OME4=0.1, OME5=0.1, OME6=0.1)}
    v_outdc2 = {'dc21': StreamGen()}
    l_outdc2 = {'dc22': StreamGen()}
    DC2_TEMP = 353
    DC2_PRESS = 760 * 2
    net1.add_component('DC2', DistillationColumn(DC2_TEMP, DC2_PRESS, DC2_recov, inletsdc2, v_outdc2, l_outdc2))
    # product outlet
    productinlet  = {'dc22': StreamGen(OME3=5, OME4=5, OME5=5, OME6=5)}
    net1.add_component('ProductOutlet', ProductRemoval(productinlet))
    # water adsorber
    ad_in = {'m31': StreamGen(MEOH=10, FA=18, H2O=0.1, OME1=0.1, OME2=0.1, OME3=0.1, OME4=0.1, OME5=0.1, OME6=0.1)}
    adwater_out = {'ad12': StreamGen()}
    adprod_out = {'ad11': StreamGen()}
    ad_recov = 1
    ad_recov_index = 3
    net1.add_component('AD1', Adsorber(ad_recov, ad_recov_index, ad_in, adprod_out, adwater_out))
    # water outlet
    waterinlet2 = {'ad12': StreamGen(H2O=20)}
    net1.add_component('WaterOutlet2', Removal(waterinlet2))
    # check stream coupling for network
    net1.check_stream_coupling()
    # test network equilibration
    net1.equilibrate_network()
