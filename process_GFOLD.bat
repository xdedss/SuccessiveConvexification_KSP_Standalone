@echo off
echo process.bat
if "%1"=="" (
echo "Empty input!"
pause
goto :end
)
@echo on
if exist %1 (rd /s/q %1)
python ./GFOLD_problem.py codegen %1
python ./postprocess.py %1 gfold_subproblem_solver
cd %1
python ./setup.py install
@echo -
@echo -
@echo OK
:end
