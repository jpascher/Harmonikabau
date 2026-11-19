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
:: Zeitstempel VOR dem if-Block setzen (innerhalb eines Blocks wird %TIMESTAMP%
:: schon beim Einlesen ersetzt und nimmt den alten Wert); PowerShell statt wmic
for /f "tokens=*" %%a in ('powershell -NoProfile -Command "Get-Date -Format 'yyyy-MM-dd HH:mm'"') do set TIMESTAMP=%%a

git diff --cached --quiet
if %errorlevel%==1 (
    git commit -m "Update vom %TIMESTAMP%"
    git push origin main
    echo Änderungen erfolgreich gepusht!
) else (
    echo Keine Änderungen zum Commit
)

pause