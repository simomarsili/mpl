rm 500.txt.scores.mtx;
../src/mpl -i 500.txt 2>/dev/null ;
if ! cmp 500.txt.scores.mtx 500.txt.SCORES.mtx >/dev/null 2>&1; then
    echo "test failed"
else
    echo "test ok"
fi
