: '
This short bash/shell script will compile and run all the ODE system examples
'

function example {
    cd ${1}
    ./run.sh
}

example examples/double-pendulum
example ../burgers
example ../lorentz
example ../star-dynamics
example ../swarmalator

cd ../..
