# Compatibility helper. Preferred: scripts/start-helixscope.ps1
$ErrorActionPreference = "Stop"
& (Join-Path $PSScriptRoot "start-helixscope.ps1") -Mode dev @args
