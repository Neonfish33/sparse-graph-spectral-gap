@echo off
cd /d "%~dp0"
echo === Toy model: exponent beta из плотности малых lambda vs alpha ===
echo (H1: alpha = 1/beta)
python toy_model.py --n 10000 --K 150 --reps 5 --out results_toy_model
echo.
echo === Cross-over: экспоненциальный хвост (theta) ===
python -c "import toy_model as T; from scipy.sparse.linalg import eigsh; print('TODO: expo sweep (theta=0.05..1.0)')"
pause
