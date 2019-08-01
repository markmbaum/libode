: '
This short bash/shell script compiles (if necessary) and runs one of the test
programs, then plots the results by running the appropriate python script.
'

#take a command line arg
test=$1
#compile test programs
make tests
#clear the output directory
python clean.py out
#run the desired test program
./bin/test_${test}.exe
#switch to scripts dir and plot the result
cd scripts
python plot_${test}.py
#back to home dir
cd ..
