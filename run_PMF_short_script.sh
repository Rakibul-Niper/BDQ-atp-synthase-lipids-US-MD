
#!/bin/bash
# ============================================
# Umbrella Sampling PMF Analysis Script
# Collects files, runs WHAM, and plots PMF
# Compatible with GROMACS 2023.5
# ============================================
export OMP_NUM_THREADS=10
# ----- User Settings -----
GMX="gmx"                  # GROMACS command
UNIT="kCal"   
TEMP="300"               # PMF output units (kJ or kCal)
BOOTSTRAP=100              # Number of bootstrap samples
BSPREFIX="bsResult"        # Bootstrap result prefix
PLOT_SCRIPT="plot_pmf.py"  # Python plotting script
OUTDIR="wham_results"      # Output directory

# ----- Make output directory -----
mkdir -p $OUTDIR
cd $OUTDIR || exit

# ----- Collect files -----
echo ">>> Collecting tpr and pullf files from ../umbrella-md/"
ls ../umbrella-md/md_umbrella*.tpr | sort -V > tpr-files.dat
ls ../umbrella-md/md_umbrella*_pullf.xvg | sort -V > pullf-files.dat

echo ">>> Number of windows found:"
wc -l tpr-files.dat pullf-files.dat

# ----- Run WHAM -----
echo ">>> Running gmx wham..."
$GMX wham -it tpr-files.dat -if pullf-files.dat \
    -o pmf.xvg -hist histo.xvg -unit $UNIT -temp ${TEMP} \
    -bsres ${BSPREFIX}.xvg -nBootstrap $BOOTSTRAP -v -b 2000 -dt 1 <<EOF
0
EOF

# Check success
if [ ! -f pmf.xvg ]; then
    echo "ERROR: pmf.xvg not generated. Check WHAM log."
    exit 1
fi

# ----- Create plotting script -----
cat > $PLOT_SCRIPT << 'EOF'
import matplotlib.pyplot as plt

# Read PMF data
x, y = [], []
with open("pmf.xvg") as f:
    for line in f:
        if line.startswith(("#", "@")):
            continue
        cols = line.split()
        if len(cols) >= 2:
            x.append(float(cols[0]))
            y.append(float(cols[1]))

# Plot
plt.figure(figsize=(7,5))
plt.plot(x, y, color="blue", linewidth=2)
plt.xlabel("Reaction Coordinate (nm)", fontsize=14)
plt.ylabel("PMF (kJ/mol)", fontsize=14)
plt.title("Potential of Mean Force", fontsize=16)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("pmf.png", dpi=300)
plt.show()
EOF

# ----- Run plot -----
echo ">>> Plotting PMF..."
python $PLOT_SCRIPT

echo ">>> Done! Results in $OUTDIR/"

