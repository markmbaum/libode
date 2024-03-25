cd build 2> /dev/null || cd ../build || mkdir -p build && cd build 2> /dev/null
cd ../build  
rm -rf ../build/* 
cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME/install .. 
cmake --build . 
make install 
cpack -G TGZ
# cpack -G RPM 
# cpack -G DEB 
# dpkg -c *.deb


