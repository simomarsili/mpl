../src/mpl -i 500.txt 2>/dev/null ;
if ! cmp 500.txt.scores 500.txt.SCORES >/dev/null 2>&1; then
    echo "test failed"
else
    echo "test ok"
fi
