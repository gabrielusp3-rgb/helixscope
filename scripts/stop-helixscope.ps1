# Stop only HelixScope processes started or recorded by start-helixscope.ps1.
# Never kills an unrelated listener merely because it occupies port 8000 or 3000.

[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$root = Split-Path -Parent $PSScriptRoot
$statePath = Join-Path $root ".helixscope-local.json"

if (-not (Test-Path $statePath)) {
    Write-Host "No .helixscope-local.json found. Nothing to stop."
    Write-Host "This script will not kill processes by port."
    exit 0
}

$state = Get-Content $statePath -Raw | ConvertFrom-Json
$stopped = @()

function Stop-RecordedPid([object]$pidValue, [bool]$reused, [string]$label) {
    if (-not $pidValue) {
        Write-Host "$label: no PID recorded."
        return
    }
    if ($reused) {
        Write-Host "$label PID $pidValue was reused (not started by the launcher). Leaving it running."
        return
    }
    try {
        $proc = Get-Process -Id ([int]$pidValue) -ErrorAction Stop
        Stop-Process -Id $proc.Id -ErrorAction Stop
        $script:stopped += "$label $($proc.Id)"
        Write-Host "Stopped $label PID $($proc.Id)"
    } catch {
        Write-Host "$label PID $pidValue is not running."
    }
}

Stop-RecordedPid $state.apiPid ([bool]$state.apiReused) "FastAPI"
Stop-RecordedPid $state.webPid ([bool]$state.webReused) "Next.js"

Remove-Item $statePath -Force -ErrorAction SilentlyContinue
if ($stopped.Count -eq 0) {
    Write-Host "No launcher-owned processes were stopped."
} else {
    Write-Host "Stopped: $($stopped -join ', ')"
}
exit 0
