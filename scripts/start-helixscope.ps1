# HelixScope local workstation launcher.
# Starts FastAPI (127.0.0.1:8000) and Next.js (http://localhost:3000).
# Does not launch Streamlit / localhost:8501.
# Does not kill unrelated processes that occupy ports 8000 or 3000.

[CmdletBinding()]
param(
    [ValidateSet("dev", "prod")]
    [string]$Mode = "prod",
    [switch]$Rebuild,
    [switch]$NoBrowser
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Write-Fail([string]$code, [string]$message) {
    Write-Host ""
    Write-Host $code
    Write-Host $message
    Write-Host "HelixScope did not start. This window can stay open so the message is visible."
    exit 1
}

function Resolve-RepoRoot {
    $here = $PSScriptRoot
    $candidate = Split-Path -Parent $here
    if (Test-Path (Join-Path $candidate "pyproject.toml")) { return $candidate }
    if (Test-Path (Join-Path $candidate "app.py") -and (Test-Path (Join-Path $candidate "apps\web\package.json"))) {
        return $candidate
    }
    throw "Could not locate the HelixScope repository root from $here"
}

function Resolve-Python([string]$root) {
    $candidates = @()
    if ($env:HELIXSCOPE_PYTHON) { $candidates += $env:HELIXSCOPE_PYTHON }
    $candidates += "C:\Python314\python.exe"
    $cmd = Get-Command python -ErrorAction SilentlyContinue
    if ($cmd) { $candidates += $cmd.Source }
    foreach ($item in $candidates) {
        if ($item -and (Test-Path $item)) { return (Resolve-Path $item).Path }
    }
    Write-Fail "DEPENDENCY MISSING" "Python interpreter not found. Set HELIXSCOPE_PYTHON to the persistent HelixScope interpreter."
}

function Test-HelixApi {
    try {
        $live = Invoke-RestMethod -Uri "http://127.0.0.1:8000/health/live" -TimeoutSec 2
        if ($live.status -eq "alive") { return "helixscope" }
        return "foreign"
    } catch {
        try {
            $null = Invoke-WebRequest -Uri "http://127.0.0.1:8000/" -UseBasicParsing -TimeoutSec 2
            return "foreign"
        } catch {
            return "empty"
        }
    }
}

function Test-HelixWeb {
    try {
        $resp = Invoke-WebRequest -Uri "http://localhost:3000/" -UseBasicParsing -TimeoutSec 3
        if ($resp.Content -match "HelixScope") { return "helixscope" }
        return "foreign"
    } catch {
        return "empty"
    }
}

function Wait-Until([scriptblock]$probe, [int]$seconds, [string]$label) {
    $deadline = (Get-Date).AddSeconds($seconds)
    while ((Get-Date) -lt $deadline) {
        if (& $probe) { return $true }
        Start-Sleep -Milliseconds 400
    }
    Write-Fail "HEALTH CHECK FAILED" "$label did not become healthy within $seconds seconds."
}

function Get-ListenPid([int]$port) {
    try {
        $row = Get-NetTCPConnection -LocalPort $port -State Listen -ErrorAction SilentlyContinue |
            Select-Object -First 1
        if ($row) { return [int]$row.OwningProcess }
    } catch { }
    return $null
}

$root = Resolve-RepoRoot
$webDir = Join-Path $root "apps\web"
$logDir = Join-Path $root ".helixscope-logs"
$statePath = Join-Path $root ".helixscope-local.json"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

if (-not (Test-Path (Join-Path $webDir "package.json"))) {
    Write-Fail "DEPENDENCY MISSING" "apps/web/package.json was not found."
}

$python = Resolve-Python $root
$node = Get-Command node -ErrorAction SilentlyContinue
$npm = Get-Command npm -ErrorAction SilentlyContinue
if (-not $node) { Write-Fail "DEPENDENCY MISSING" "Node.js is not on PATH." }
if (-not $npm) { Write-Fail "DEPENDENCY MISSING" "npm is not on PATH." }

if (-not (Test-Path (Join-Path $webDir "node_modules"))) {
    Write-Fail "DEPENDENCY MISSING" "apps/web/node_modules is missing. Run: cd apps\web; npm ci"
}

$env:PYTHONPATH = "$root;$root\services\api"
$envFile = Join-Path $webDir ".env.local"
if (-not (Test-Path $envFile)) {
    Copy-Item (Join-Path $webDir ".env.example") $envFile
}

$apiState = Test-HelixApi
$webState = Test-HelixWeb
$apiPid = $null
$webPid = $null
$apiReused = $false
$webReused = $false

if ($apiState -eq "foreign") {
    Write-Fail "PORT OCCUPIED" "Port 8000 is in use by a process that is not HelixScope FastAPI (/health/live). Refusing to kill it."
}
if ($webState -eq "foreign") {
    Write-Fail "PORT OCCUPIED" "Port 3000 is in use by a process that is not HelixScope Next.js. Refusing to kill it."
}

if ($apiState -eq "helixscope") {
    $apiReused = $true
    $apiPid = Get-ListenPid 8000
    Write-Host "Reusing HelixScope FastAPI on 127.0.0.1:8000 (PID $($apiPid))"
} else {
    $apiArgs = @("-m", "uvicorn", "helixscope_api.main:app", "--host", "127.0.0.1", "--port", "8000")
    if ($Mode -eq "dev") { $apiArgs += "--reload" }
    $apiOut = Join-Path $logDir "fastapi.out.log"
    $apiErr = Join-Path $logDir "fastapi.err.log"
    try {
        $apiProc = Start-Process -FilePath $python -ArgumentList $apiArgs -WorkingDirectory $root -PassThru -WindowStyle Hidden -RedirectStandardOutput $apiOut -RedirectStandardError $apiErr
    } catch {
        Write-Fail "FASTAPI FAILED" $_.Exception.Message
    }
    $apiPid = $apiProc.Id
    Wait-Until { (Test-HelixApi) -eq "helixscope" } 60 "FastAPI /health/live" | Out-Null
    try {
        $null = Invoke-RestMethod -Uri "http://127.0.0.1:8000/health/ready" -TimeoutSec 5
    } catch {
        Write-Fail "HEALTH CHECK FAILED" "FastAPI /health/ready did not succeed. See $apiErr"
    }
    Write-Host "FastAPI started (PID $apiPid)"
}

if ($webState -eq "helixscope") {
    $webReused = $true
    $webPid = Get-ListenPid 3000
    Write-Host "Reusing HelixScope Next.js on http://localhost:3000 (PID $($webPid))"
} else {
    $nextDir = Join-Path $webDir ".next"
    if ($Mode -eq "prod") {
        $needBuild = $Rebuild -or -not (Test-Path $nextDir)
        if ($needBuild) {
            Write-Host "Building Next.js production bundle..."
            Push-Location $webDir
            try {
                & npm run build
                if ($LASTEXITCODE -ne 0) { Write-Fail "NEXT.JS FAILED" "npm run build failed." }
            } finally {
                Pop-Location
            }
        } else {
            Write-Host "Using existing apps/web/.next (pass -Rebuild to force a new production build)."
        }
    }
    $webOut = Join-Path $logDir "next.out.log"
    $webErr = Join-Path $logDir "next.err.log"
    $nextBin = Join-Path $webDir "node_modules\next\dist\bin\next"
    if (-not (Test-Path $nextBin)) {
        Write-Fail "DEPENDENCY MISSING" "Next.js binary missing at $nextBin. Run npm ci in apps/web."
    }
    $nodeExe = $node.Source
    $nextArgs = @($nextBin)
    if ($Mode -eq "prod") {
        $nextArgs += @("start", "--port", "3000", "--hostname", "localhost")
    } else {
        $nextArgs += @("dev", "--port", "3000", "--hostname", "localhost")
    }
    try {
        $webProc = Start-Process -FilePath $nodeExe -ArgumentList $nextArgs -WorkingDirectory $webDir -PassThru -WindowStyle Hidden -RedirectStandardOutput $webOut -RedirectStandardError $webErr
    } catch {
        Write-Fail "NEXT.JS FAILED" $_.Exception.Message
    }
    $webPid = $webProc.Id
    Wait-Until { (Test-HelixWeb) -eq "helixscope" } 90 "Next.js http://localhost:3000" | Out-Null
    Write-Host "Next.js started (PID $webPid) mode=$Mode"
}

try {
    $ready = Invoke-RestMethod -Uri "http://127.0.0.1:8000/health/ready" -TimeoutSec 5
    if ($ready.status -ne "ready") {
        Write-Fail "HEALTH CHECK FAILED" "FastAPI /health/ready status is '$($ready.status)'."
    }
} catch {
    Write-Fail "HEALTH CHECK FAILED" "Frontend host cannot reach FastAPI /health/ready."
}

$identity = $null
try {
    $identity = Invoke-RestMethod -Uri "http://127.0.0.1:8000/api/v1/system/identity" -TimeoutSec 5
} catch {
    Write-Fail "HEALTH CHECK FAILED" "FastAPI /api/v1/system/identity is not readable. Restart HelixScope with a current API process."
}
$contractPath = Join-Path $webDir "lib\api\generated\contract.ts"
if (-not (Test-Path $contractPath)) {
    Write-Fail "DEPENDENCY MISSING" "apps/web/lib/api/generated/contract.ts is missing. Run npm run api:types in apps/web."
}
$contractText = Get-Content $contractPath -Raw
$expectedHash = $null
if ($contractText -match 'EXPECTED_OPENAPI_SHA256 = "([a-f0-9]{64})"') {
    $expectedHash = $Matches[1]
}
if (-not $expectedHash) {
    Write-Fail "DEPENDENCY MISSING" "Could not read EXPECTED_OPENAPI_SHA256 from contract.ts."
}
if ($identity.openapi_sha256 -ne $expectedHash) {
    Write-Fail "API/FRONTEND MISMATCH" "Restart HelixScope using the supported launcher after regenerating API types. This is a software compatibility state, not a scientific result. expected=$expectedHash running=$($identity.openapi_sha256)"
}
Write-Host "API identity ok product=$($identity.product_version) core=$($identity.core_version) hash=$($identity.openapi_sha256.Substring(0,12))"

$state = [ordered]@{
    startedAt = (Get-Date).ToString("o")
    mode = $Mode
    apiPid = $apiPid
    webPid = $webPid
    apiReused = $apiReused
    webReused = $webReused
    python = $python
    fastapi = "http://127.0.0.1:8000"
    next = "http://localhost:3000"
    docs = "http://127.0.0.1:8000/docs"
    streamlitLegacy = "http://localhost:8501 (not started)"
    productVersion = $identity.product_version
    coreVersion = $identity.core_version
    openapiSha256 = $identity.openapi_sha256
}
$state | ConvertTo-Json | Set-Content -Path $statePath -Encoding utf8

Write-Host ""
Write-Host "PRIMARY PRODUCT PLATFORM: NEXT.JS + FASTAPI + HELIXSCOPE CORE"
Write-Host "FastAPI:  http://127.0.0.1:8000"
Write-Host "Ready:    http://127.0.0.1:8000/health/ready"
Write-Host "Docs:     http://127.0.0.1:8000/docs"
Write-Host "Next.js:  http://localhost:3000"
Write-Host "FastAPI PID: $apiPid  (reused=$apiReused)"
Write-Host "Next.js PID: $webPid  (reused=$webReused)"
Write-Host "State file: $statePath"
Write-Host "Stop with:  powershell -File scripts\stop-helixscope.ps1"
Write-Host "Streamlit was not started."

if (-not $NoBrowser) {
    Start-Process "http://localhost:3000"
}
exit 0
