# python fitting_model_params.py --N 100 -df 0 --fitting &
# python fitting_model_params.py --N 100 -df 1 --fitting &
# python fitting_model_params.py --N 100 -df 2 --fitting &
# python fitting_model_params.py --N 100 -df 3 --fitting &
# python fitting_model_params.py --N 100 -df 4 --fitting &
# python fitting_model_params.py --N 100 -df 5 --fitting &
# python fitting_model_params.py --N 100 -df 6 --fitting &
# python fitting_model_params.py --N 100 -df 7 --fitting &
# python fitting_model_params.py --N 100 -df 8 --fitting &
# python fitting_model_params.py --N 100 -df 9 --fitting &
# python fitting_model_params.py --N 100 -df 10 --fitting &
# python fitting_model_params.py --N 100 -df 11 --fitting &
# python fitting_model_params.py --N 100 -df 12 --fitting &

# ########################################################
for i in {0..13}
do
    echo sget '../ring_model/data/fitting_data_'$i'.npy'
done
echo ''
