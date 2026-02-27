@echo off
chcp 65001
cd /d %~dp0

echo === Harmonikabau Repository Update ===

:: .gitignore prüfen (falls nicht vorhanden)
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
    git commit -m "Update vom %date% %time%"
    git push origin master
    echo Änderungen erfolgreich gepusht!
) else (
    echo Keine Änderungen zum Commit
)

pause