@echo off
chcp 65001
cd /d %~dp0

echo === Harmonikabau Repository Update ===

:: .gitignore prüfen
if not exist .gitignore (
    echo *.aux > .gitignore
    echo *.log >> .gitignore
    echo *.out >> .gitignore
    echo *.synctex.gz >> .gitignore
    echo .gitignore erstellt
)

:: Status anzeigen
git status

:: Änderungen hinzufügen
git add .

:: Prüfen ob Änderungen vorhanden
git diff --cached --quiet
if %errorlevel%==1 (
    for /f "tokens=*" %%a in ('wmic os get localdatetime ^| find "."') do set datetime=%%a
    set TIMESTAMP=%datetime:~0,4%-%datetime:~4,2%-%datetime:~6,2% %datetime:~8,2%:%datetime:~10,2%
    git commit -m "Update vom %TIMESTAMP%"
    git push origin main
    echo Änderungen erfolgreich gepusht!
) else (
    echo Keine Änderungen zum Commit
)

pause