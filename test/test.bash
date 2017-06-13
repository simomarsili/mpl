rm test.txt.scores.mtx;
../src/mpl -i test.txt 2>/dev/null ;
if ! cmp test.txt.scores.mtx test.txt.SCORES.mtx >/dev/null 2>&1; then
    echo "test failed"
else
    echo "test ok"
fi
