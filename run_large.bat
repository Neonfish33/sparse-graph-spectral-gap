@echo off
cd /d "%~dp0"
echo === Large-n test: d_min=2 vs 3, n=3e4 and 1e5 ===
echo (можно запускать повторно: --resume продолжит с места обрыва)
python run_experiments.py --n 30000 100000 --gammas 2.1 2.5 3.0 --dmin 2 3 ^
  --no-planted --reps 100 --dmax 50 --jobs 6 --save-curves --resume --out results_large
echo === Analysis ===
python analyze_results.py --csv results_large\diagnostics.csv --out results_large\analysis
echo === DONE. See results_large\analysis\core_degree_comparison.csv and summary.txt ===
pause
