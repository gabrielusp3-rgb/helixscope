"use client";

export default function ModuleError({
  error,
  reset,
}: {
  error: Error & { digest?: string };
  reset: () => void;
}) {
  return (
    <div className="hs-page">
      <h1>Module failed</h1>
      <p className="error-box">{error.message}</p>
      <p className="hs-lede">The workstation shell remains available.</p>
      <button type="button" className="btn-primary" onClick={reset}>
        Retry module
      </button>
    </div>
  );
}
