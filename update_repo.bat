@echo off
chcp 65001
cd /d %~dp0

:: Überprüfe, ob .gitignore existiert, und erstelle eine, falls nicht vorhanden
if not exist .gitignore (
    echo *.aux > .gitignore
    echo *.log >> .gitignore
    echo *.out >> .gitignore
    echo *.synctex.gz >> .gitignore
    echo .gitignore erstellt mit Standard-LaTeX-Ausschlüssen.
)

:: Überprüfe Git-Status
git status

:: Stage alle Änderungen im gesamten Repository
git add .

:: Stelle sicher, dass das Unterverzeichnis explizit geaddet wird (inklusive neuer Dateien)


:: Prüfe, ob es gestagte Änderungen gibt
git diff --cached --quiet
if %errorlevel%==1 (
    :: Commit mit detaillierter Nachricht
    git commit -m "Automatischer Commit: Änderungen im gesamten Repository und explizites Update für 2/html/ vom %date% %time%"
    :: Push mit Fehlerbehandlung
    git push origin main || (
        echo Fehler beim Pushen zum Repository.
        pause
        exit /b 1
    )
    echo Änderungen erfolgreich committet und gepusht.
) else (
    echo Keine Änderungen zum Commit vorhanden.
)

pause