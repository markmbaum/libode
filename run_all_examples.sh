: '
This short bash/shell script will compile and run all the ODE system examples
'

function example {
    cwd=$(pwd)
    cd ${1}
    ./run.sh
    cd $cwd
}

example examples/burgers
example examples/double-pendulum
example examples/lorentz
example examples/nonlinear-diffusion
example examples/star-dynamics
example examples/swarmalator
