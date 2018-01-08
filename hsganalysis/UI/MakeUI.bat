@ECHO OFF
set input="%~1"
set postpend="_ui.py"
set output=%input:~1,-4%
pyuic5 %input% -o "%input:~0,-4%_ui.py"