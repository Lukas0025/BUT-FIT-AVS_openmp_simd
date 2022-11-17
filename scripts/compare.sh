

SCRIPT_ROOT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
CALCULATORS=("ref" "line" "batch")
VALID=1

for size in {1..4096}; do
    for calc in "${CALCULATORS[@]}"; do
        ./mandelbrot -s $size -i 100 -c $calc --batch cmp_$calc.npz
    done

    echo "Reference vs line $size"
    python3 ${SCRIPT_ROOT_PATH}/compare.py cmp_ref.npz cmp_line.npz  || VALID=0

    echo "Reference vs batch $size"
    python3 ${SCRIPT_ROOT_PATH}/compare.py cmp_ref.npz cmp_batch.npz || VALID=0


    echo "Batch vs line $size"
    python3 ${SCRIPT_ROOT_PATH}/compare.py cmp_line.npz cmp_batch.npz || VALID=0

    if [ "$VALID" -eq 0 ]; then
        echo "Test fail $size";
        exit 1
    fi
done

if [ "$VALID" -eq 1 ]; then
    echo "Test passed";
else
    echo "Test failed";
    exit 1
fi



