# HelixScope web (Prompt 4)

Next.js App Router workstation. It does not compute biology. HelixScope API v1 is the scientific authority.

Internal package version `0.4.0-prompt4` is not a HelixScope product release. Product version remains **0.24.3-19**.

Streamlit (`app.py`) is a retired reference UI. Run this workstation plus FastAPI.

## Requirements

- Node.js 20.9+ (this machine used Node 24 LTS)
- npm (do not add yarn.lock or pnpm-lock.yaml)
- FastAPI from Prompt 2 on `http://127.0.0.1:8000`
- Browser origin `http://localhost:3000` (FastAPI CORS default)

## Interpreter / API

Keep the persistent HelixScope Python interpreter for Core/API work (`C:\Python314\python.exe` on the development machine). Do not select disposable venvs such as `%TEMP%\hs-core-clean` or `.venv-api-prompt2`.

## Run

Preferred (repository root):

```powershell
powershell -File ..\..\scripts\start-helixscope.ps1
```

That starts FastAPI and Next.js together and opens http://localhost:3000.

Manual fallback:

Terminal 1 (repository root):

```powershell
$env:PYTHONPATH = "$((Get-Location).Path);$((Get-Location).Path)\services\api"
python -m uvicorn helixscope_api.main:app --host 127.0.0.1 --port 8000
```

Terminal 2:

```powershell
cd apps/web
copy .env.example .env.local
npm ci
npm run api:types
npm run mediapipe:assets
npm run dev
```

Production preview:

```powershell
npm run build
npm run start
```

## Scripts

- `npm run api:types` — generate `lib/api/generated/` from `services/api/openapi.json`
- `npm run api:types:check` — fail if OpenAPI changed without regeneration
- `npm run typecheck`
- `npm run lint`
- `npm run test` — Vitest
- `npm run test:e2e` — Playwright against `next start` + live FastAPI
- `npm run mediapipe:assets` — copy WASM; fetch pinned Hand Landmarker `/1/` model

## Contract

Never edit `lib/api/generated/`. After FastAPI OpenAPI changes, regenerate types.
