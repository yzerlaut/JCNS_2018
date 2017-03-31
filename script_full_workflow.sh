
##########################################################################
##########  The paper configuration                        ###############
##########################################################################
# we start with transfer functions
cd transfer_functions
# FS cell
python tf_simulation.py FS-cell CONFIG1 -s --SEED 4 --tstop 10
python theoretical_tools.py -f data/FS-cell_CONFIG1.npy --With_Square
# RS cell
python tf_simulation.py RS-cell CONFIG1 -s --SEED 4 --tstop 10
python theoretical_tools.py -f data/RS-cell_CONFIG1.npy --With_Square
cd ..
# now network simulations
cd network_simulations
python ntwk_sim_demo.py --CONFIG RS-cell--FS-cell--CONFIG1 -f data/config1.py
# now mean field prediction, passing the output of ntwk sim for comparison
python compare_with_mean_field.py RS-cell--FS-cell--CONFIG1 -f data/config1.py
cd ..

# ##########################################################################
# ##########  A new configuration                        ###############
# ##########################################################################
# # we start with transfer functions
# cd transfer_functions
# # FS cell
# python tf_simulation.py FS-cell2 CONFIG1 -s --SEED 4 --tstop 10
# python theoretical_tools.py -f data/FS-cell2_CONFIG1.npy --With_Square
# # RS cell2
# python tf_simulation.py RS-cell2 CONFIG1 -s --SEED 4 --tstop 10
# python theoretical_tools.py -f data/RS-cell2_CONFIG1.npy --With_Square
# cd ..
# # now network simulations
# cd network_simulations
# python ntwk_sim_demo.py --CONFIG RS-cell2--FS-cell2--CONFIG1 -f data/config2.py
# # now mean field prediction, passing the output of ntwk sim for comparison
# python compare_with_mean_field.py RS-cell2--FS-cell2--CONFIG1 -f data/config2.py
# cd ..
