: '
This short bash/shell script runs all the tests!
'

function test {

    printf "\n--------------------------------------------------------------------------------\n"
    echo "RUNNING TEST:" ${1}
    printf "\n"
    ./bin/test_${1}.exe
}

for test in adapt all conv fixed linalg newton snaps stiff work
do
    test $test
done
