# #!/bin/bash



# RUN_DIR="/sps/nemo/snemo/snemo_data/reco_data/UDD/delta-tdc-10us-v3/"
# PROCESSED_FILE="/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durations_betabeta.txt"



# for RUN in $(seq 2204 2570)
# do
#     if ! awk '{print $1}' "$PROCESSED_FILE" | grep -qx "$RUN"; then
#         if [ -f "$RUN_DIR/snemo_run-$RUN_udd.root" ]; then
#             echo "New run detected : $RUN"
#             root -l -b -q "/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durationner.cpp($RUN)"
#         fi
#     fi
# done

#!/bin/bash
RUN_DIR="/sps/nemo/snemo/snemo_data/reco_data/UDD/delta-tdc-10us-v3/"
PROCESSED_FILE="/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durations_betabeta.txt"
for RUN in $(seq 2204 2570)
do
    if ! awk '{print $1}' "$PROCESSED_FILE" | grep -qx "$RUN"; then
        echo "New run detected : $RUN"
        root -l -b -q "/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durationner.cpp($RUN)"
    fi
done