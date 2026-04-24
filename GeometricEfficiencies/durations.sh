#!/bin/bash
# -*- coding: utf-8 -*-
# calc_run_duration.sh
#
# Розраховує сумарну кількість секунд активного датасету
# по номерам ранів від 2204 до 2570.
# Якщо номер рану відсутній у файлі — вважається тривалістю 1800 секунд.
#
# Використання:
#   bash calc_run_duration.sh run_durations_betabeta.txt

PROCESSED_FILE="${1:-/sps/nemo/scratch/kfilonen/Calibration/SourcetoOM/Timer/run_durations_betabeta.txt}"

RUN_START=2204
RUN_END=2570
DEFAULT_DURATION=1800

if [ ! -f "$PROCESSED_FILE" ]; then
    echo "ERROR: File not found: $PROCESSED_FILE"
    exit 1
fi

total=0
found=0
missing=0

for RUN in $(seq $RUN_START $RUN_END); do
    # Шукаємо рядок з цим номером рану та трьома колонками (RUN STARTTIME DURATION)
    DURATION=$(awk -v run="$RUN" '$1 == run && NF >= 3 { print $3 }' "$PROCESSED_FILE")

    if [ -n "$DURATION" ]; then
        total=$((total + DURATION))
        found=$((found + 1))
    else
        total=$((total + DEFAULT_DURATION))
        missing=$((missing + 1))
    fi
done

total_runs=$((RUN_END - RUN_START + 1))
hours=$(echo "scale=2; $total / 3600" | bc)
days=$(echo "scale=2; $total / 86400" | bc)

echo "Діапазон ранів     : $RUN_START – $RUN_END"
echo "Всього ранів       : $total_runs"
echo "  знайдено у файлі : $found"
echo "  відсутніх (×${DEFAULT_DURATION}с): $missing"
echo "─────────────────────────────────────"
echo "Сумарна тривалість : $total секунд"
echo "                   ≈ $hours годин"
echo "                   ≈ $days діб"
