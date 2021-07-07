for th in 0 60 120 180 240 270
do
    mkdir th$th
    cd th$th
    cp ../Ffmulti.inp .
    sed -i s/pth/$th/g Ffmulti.inp
    ../Fracflow Ffmulti
    cd ..
done
