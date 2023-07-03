for i in 16 32 64 128
do
	for j in 16 32 64 128
	do
		mkdir "Np$i"_"Nt$j"
		cd "Np$i"_"Nt$j"
		../a.out $i $j config.tether
		mpirun -np 4 python3 ../batchrun.py ../in.setup
		cd ..
	done
done

# Sets up the directories and such for
