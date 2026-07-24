$ErrorActionPreference = "Stop"
Push-Location frontend
try {
    npm install
    npm run typecheck
    npm run build
}
finally {
    Pop-Location
}
Write-Host "Frontend build completed. Run: python app.py"
