cd ../transfer_functions/
python tf_simulation.py FS-cell CONFIG1 -s
python theoretical_tools.py -f data/FS-cell_CONFIG1.npy --With_Square
python tf_simulation.py RS-cell CONFIG1 -s
python theoretical_tools.py -f data/RS-cell_CONFIG1.npy --With_Square
cd ../ntwk_analysis/
python phase_space.py RS-cell FS-cell CONFIG1 --ext_drive 0.
