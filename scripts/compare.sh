

SCRIPT_ROOT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
CALCULATORS=("ref" "line" "batch")
SIZES=("32" "64" "128" "92" "93" "1039" "256" "1024" "4096")
VALID=1

for size in "${SIZES[@]}"; do
    for calc in "${CALCULATORS[@]}"; do
        ./mandelbrot -s $size -i 100 -c $calc --batch cmp_$calc.npz
    done

    echo "Reference vs line $size"
    python3 ${SCRIPT_ROOT_PATH}/compare.py cmp_ref.npz cmp_line.npz  || VALID=0

    echo "Reference vs batch $size"
    python3 ${SCRIPT_ROOT_PATH}/compare.py cmp_ref.npz cmp_batch.npz || VALID=0


    echo "Batch vs line $size"
    python3 ${SCRIPT_ROOT_PATH}/compare.py cmp_line.npz cmp_batch.npz || VALID=0
done

if [ "$VALID" -eq 1 ]; then
    echo "Test passed";
else
    echo "Test failed";
    exit 1
fi



