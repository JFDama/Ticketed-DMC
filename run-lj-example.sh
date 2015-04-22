# Check for the expected input parameters.
if [ $# -ne 3 ]; then
	echo "Usage: $0 gamma lambda n_replicates"
	echo "All parameters are required."
# Run the example if three are present.
else
	# Create a temporary file.
	cp lj-example.cpp lj-ex-tmp.cpp
	# Change the parameters of the temporary.
	perl -pi -e 's/paper_gamma = [0-9.]+/paper_gamma = '"$1"'/' lj-ex-tmp.cpp
	perl -pi -e 's/paper_lambda = [0-9.]+/paper_lambda = '"$2"'/' lj-ex-tmp.cpp
	perl -pi -e 's/n_runs = [0-9]+/n_runs = '"$3"'/' lj-ex-tmp.cpp
	# Compile the temporary and then remove it.
	clang++ --std=c++11 -O3 lj-ex-tmp.cpp -o ljdmc.x
	# Run the compiled executable and then remove it.
	outfile=${1}-${2}-${3}-out.txt
	./ljdmc.x > $outfile
	rm -f ljdmc.x
	rm -f lj-ex-tmp.cpp
fi
